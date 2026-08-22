function [coeffs, info] = filtercoeff(spec)
% FILTERCOEFF  Faithful MATLAB port of the FilterCoeff C99 library.
%
%   [coeffs, info] = filtercoeff(spec)
%
% Turns a filter specification into ready-to-use coefficients. This is a
% direct port of the C algorithms in src/fce_math.c, src/fce_fir.c and
% src/fce_iir.c, so it produces the same coefficients as the C library and
% the Python port (python/filtercoeff.py) to ~1e-15.
%
%   spec.kind   : 'fir' | 'iir'
%   spec.fs, spec.fc1, spec.fc2 : sample rate and cutoffs / band edges [Hz]
%   FIR:
%     spec.fir_type : 'lowpass'|'highpass'|'bandpass'|'bandstop'
%                     |'hilbert'|'differentiator'
%     spec.num_taps : tap count (0 = auto via Kaiser)
%     spec.window   : 'rectangular'|'hann'|'hamming'|'blackman'|'kaiser'
%                     |'blackman_harris'|'bartlett'|'tukey'
%     spec.kaiser_beta : 0 = auto from stopband_atten_db
%     spec.normalization : 'auto'|'dc'|'nyquist'|'passband_peak'|'none'
%   IIR:
%     spec.iir_family : 'butterworth'|'chebyshev1'|'chebyshev2'
%                       |'elliptic'|'bessel'
%     spec.iir_type   : 'lowpass'|'highpass'|'bandpass'|'bandstop'
%     spec.order      : prototype order (0 = auto)
%     spec.passband_ripple_db, spec.stopband_atten_db : dB
%     spec.edge1_hz, spec.edge2_hz : opposite-band edges (auto-order only)
%     spec.sos_order : 'default'|'pole_radius_asc'|'pole_radius_desc'
%                      |'internal_gain'
%
%   coeffs : FIR -> 1xN coefficient vector
%            IIR -> 1x(5*S) flattened SOS, layout per section
%                   {b0, b1, b2, a1, a2}  (a0 == 1 implicit)
%   info   : struct with metadata (num_taps / order / num_sections /
%            design_fc1 / design_fc2 / kaiser_beta / ...).
%
% Sign convention (IIR, SOS):  a0 == 1 implicit, a1/a2 negated in the
% recurrence y[n] = b0 x[n] + b1 x[n-1] + b2 x[n-2] - a1 y[n-1] - a2 y[n-2].
%
% See also FILTERCOEFF_QUANTIZE, FILTERCOEFF_DEMO, VALIDATE_FILTERCOEFF.

    s = filtercoeff_defaults(spec);

    if strcmp(s.kind, 'fir')
        n = s.num_taps;
        if n == 0
            if ~strcmp(s.window, 'kaiser') || s.stopband_atten_db <= 0 ...
                    || s.transition_hz <= 0
                error('FilterCoeff:spec', ...
                      'auto taps require Kaiser + attenuation + transition');
            end
            n = fce_kaiser_taps(s.stopband_atten_db, s.transition_hz, s.fs);
            n = min(n, 2048);   % FCE_MAX_FIR_TAPS
        end
        s.num_taps = n;
        r = fce_fir_design(s);
        coeffs = r.h;
        info = r;
        info.kind = 'fir';
    else
        order = s.order;
        if order == 0
            a = fce_iir_auto_order(s);
            order = a.order;
            s.design_fc1 = a.design_fc1;
            s.design_fc2 = a.design_fc2;
        else
            if ~isfield(s, 'design_fc1') || isempty(s.design_fc1)
                s.design_fc1 = s.fc1;
            end
            if (~isfield(s, 'design_fc2') || isempty(s.design_fc2)) ...
                    && isfield(s, 'fc2')
                s.design_fc2 = s.fc2;
            end
        end
        s.order = order;
        r = fce_iir_design(s);
        coeffs = zeros(1, 5 * r.num_sections);
        for k = 1:r.num_sections
            coeffs((k-1)*5 + (1:5)) = r.sos(k, :);
        end
        info = r;
        info.kind = 'iir';
    end
end

% ======================================================================
% defaults
% ======================================================================
function s = filtercoeff_defaults(spec)
    s = spec;
    if ~isfield(s, 'fs') || s.fs <= 0
        error('FilterCoeff:spec', 'fs must be > 0');
    end
    if ~isfield(s, 'fc1'), s.fc1 = 0; end
    if ~isfield(s, 'fc2'), s.fc2 = 0; end
    if ~isfield(s, 'normalization'), s.normalization = 'auto'; end
    if ~isfield(s, 'precision'), s.precision = 'float64'; end
    if ~isfield(s, 'window'), s.window = 'hamming'; end
    if ~isfield(s, 'num_taps'), s.num_taps = 0; end
    if ~isfield(s, 'kaiser_beta'), s.kaiser_beta = 0; end
    if ~isfield(s, 'order'), s.order = 0; end
    if ~isfield(s, 'passband_ripple_db'), s.passband_ripple_db = 0; end
    if ~isfield(s, 'stopband_atten_db'), s.stopband_atten_db = 0; end
    if ~isfield(s, 'transition_hz'), s.transition_hz = 0; end
    if ~isfield(s, 'edge1_hz'), s.edge1_hz = 0; end
    if ~isfield(s, 'edge2_hz'), s.edge2_hz = 0; end
    if ~isfield(s, 'sos_order'), s.sos_order = 'default'; end
    if ~isfield(s, 'tukey_alpha'), s.tukey_alpha = 0.5; end
end

% ======================================================================
% math core (port of fce_math.c)
% ======================================================================
function y = fce_sinc(x)
    if x == 0
        y = 1;
    else
        px = pi * x;
        y = sin(px) / px;
    end
end

function y = fce_i0(x)
    s = 1; term = 1; x2 = 0.25 * x * x; k = 1;
    while k < 200
        term = term * x2 / (k * k);
        s = s + term;
        if term < 1e-18 * s
            break;
        end
        k = k + 1;
    end
    y = s;
end

function y = fce_ellipk(m)
    if m <= 0
        y = pi / 2;
    elseif m >= 1
        y = 1e300;
    else
        a = 1; b = sqrt(1 - m);
        while true
            t = 0.5 * (a + b);
            b = sqrt(a * b);
            a = t;
            if a - b < 1e-15 * a
                break;
            end
        end
        y = pi / (2 * a);
    end
end

function y = fce_ellipkm1(m)
    if m <= 0
        y = 1e300;
    elseif m >= 1
        y = pi / 2;
    else
        a = 1; b = sqrt(m);
        while true
            t = 0.5 * (a + b);
            b = sqrt(a * b);
            a = t;
            if a - b < 1e-15 * a
                break;
            end
        end
        y = pi / (2 * a);
    end
end

function [sn, cn, dn] = fce_ellipj(u, m)
    if m <= 0
        sn = sin(u); cn = cos(u); dn = 1;
        return;
    end
    if m >= 1
        sn = tanh(u); cn = 1 / cosh(u); dn = 1 / cosh(u);
        return;
    end
    a = 1; b = sqrt(1 - m); c = sqrt(m);
    n = 0;
    A = zeros(1, 64); C = zeros(1, 64);
    A(1) = a; C(1) = c;
    while n < 63
        an = 0.5 * (A(n+1) + b);
        cn = 0.5 * (A(n+1) - b);
        b = sqrt(A(n+1) * b);
        A(n+2) = an;
        C(n+2) = cn;
        n = n + 1;
        if C(n+1) < 1e-16 * A(n+1)
            break;
        end
    end
    phi = A(n+1) * 2^n * u;
    for i = n:-1:1
        t = (C(i+1) / A(i+1)) * sin(phi);
        if t > 1, t = 1; end
        if t < -1, t = -1; end
        phi = 0.5 * (asin(t) + phi);
    end
    sn = sin(phi);
    cn = cos(phi);
    dn = sqrt(1 - m * sn * sn);
end

function kp = fce_cx_complement(kx)
    kp = sqrt((1 - kx) * (1 + kx));
end

function u = fce_arc_jac_sn(w, m)
    k = sqrt(m);
    if k > 1
        u = nan(1, 1) + 1i * nan(1, 1);
        return;
    end
    if k == 1
        u = atanh(w);
        return;
    end
    ks = k;
    niter = 0;
    while ks(niter+1) ~= 0
        kp = fce_cx_complement(ks(niter+1));
        ks(niter+2) = (1 - kp) / (1 + kp);
        niter = niter + 1;
        if niter >= 62
            u = nan(1,1) + 1i * nan(1,1);
            return;
        end
    end
    K = 1;
    for i = 1:niter
        K = K * (1 + ks(i+1));
    end
    K = K * pi / 2;
    wns = w;
    for i = 0:niter-1
        num = 2 * wns(i+1);
        arg = wns(i+1) * ks(i+1);
        cc = sqrt(1 - arg * arg);
        den = (1 + ks(i+2)) * (1 + cc);
        wns(i+2) = num / den;
        if ~(isfinite(real(wns(i+2))) && isfinite(imag(wns(i+2))))
            u = nan(1,1) + 1i * nan(1,1);
            return;
        end
    end
    u = (2 / pi) * asin(wns(niter+1)) * K;
end

function y = fce_ellipdeg(n, m1)
    K1 = fce_ellipk(m1);
    K1p = fce_ellipkm1(m1);
    q1 = exp(-pi * K1p / K1);
    q = q1 ^ (1 / n);
    num = 0;
    for k = 0:8
        num = num + q ^ (k * (k + 1));
    end
    den = 1;
    for k = 1:9
        den = den + 2 * q ^ (k * k);
    end
    y = 16 * q * (num / den) ^ 4;
end

function y = fce_polyval(x, c)
    y = 0;
    for i = numel(c):-1:1
        y = y * x + c(i);
    end
end

function r = fce_roots_durand_kerner(coeff, deg, start)
    lead = coeff(deg + 1);
    if lead == 0
        error('FilterCoeff:numerical', 'leading coefficient zero');
    end
    a = coeff / lead;   % 1 x (deg+1)
    if nargin >= 3 && ~isempty(start)
        r = start;
    else
        sc = sum(abs(a)) / (deg + 1);
        if sc <= 0
            error('FilterCoeff:numerical', 'zero polynomial');
        end
        r = zeros(1, deg);
        for i = 1:deg
            ang = 2 * pi * (i - 1) / deg;
            rad = 0.7 + 0.3 * mod((i - 1) * 7, 5) / 4.0;
            r(i) = rad * exp(1i * ang) * sc;
        end
    end
    converged = false;
    for it = 1:200
        converged = true;
        for i = 1:deg
            pv = 1 + 0i; dp = 0 + 0i;
            for j = deg:-1:1
                dp = dp * r(i) + pv;
                pv = pv * r(i) + a(j);
            end
            ssum = 0 + 0i;
            for j = 1:deg
                if j ~= i
                    ssum = ssum + 1 / (r(i) - r(j));
                end
            end
            w = pv / dp;
            den = 1 - w * ssum;
            if abs(den) < 1e-300
                error('FilterCoeff:numerical', 'denominator too small');
            end
            step = w / den;
            r(i) = r(i) - step;
            if abs(step) > 1e-12 * (1 + abs(r(i)))
                converged = false;
            end
        end
        if converged
            break;
        end
    end
    if ~converged
        error('FilterCoeff:numerical', 'root finder did not converge');
    end
end

function z = fce_campos_zeros(n)
    if n == 1
        z = -1;
        return;
    end
    s = fce_polyval(n, [0, 0, 2, 0, -3, 1]);
    b3 = fce_polyval(n, [16, -8]) / s;
    b2 = fce_polyval(n, [-24, -12, 12]) / s;
    b1 = fce_polyval(n, [8, 24, -12, -2]) / s;
    b0 = fce_polyval(n, [0, -6, 0, 5, -1]) / s;
    r = fce_polyval(n, [0, 0, 2, 1]);
    a1 = fce_polyval(n, [-6, -6]) / r;
    a2 = 6 / r;
    z = zeros(1, n);
    for k = 1:n
        kk = k;
        z(k) = (a2 * kk * kk + a1 * kk) ...
               + 1i * (((b3 * kk + b2) * kk + b1) * kk + b0);
    end
end

function poles = fce_bessel_poles(n)
    a = zeros(1, n + 1);
    for k = 0:n
        num = 1;
        for i = 0:(2*n - k - 1)
            num = num * (2 * n - k - i);
        end
        den = 2 ^ (n - k);
        for i = 2:k
            den = den * i;
        end
        for i = 2:(n - k)
            den = den * i;
        end
        a(k + 1) = num / den;
    end
    yz = fce_campos_zeros(n);
    start = 1 ./ yz;
    poles = fce_roots_durand_kerner(a, n, start);
    ssum = sum(real(poles));
    if abs(ssum + a(n) / a(n + 1)) > 1e-8 * (1 + abs(ssum))
        error('FilterCoeff:numerical', 'bessel poles sanity failed');
    end
end

% ======================================================================
% biquad helpers
% ======================================================================
function h = fce_eval_biquad(c, w)
    % c = [b0, b1, b2, a1, a2]
    cw = cos(w); sw = sin(w);
    c2w = cos(2 * w); s2w = sin(2 * w);
    num = (c(1) + c(2) * cw + c(3) * c2w) ...
          + 1i * (-(c(2) * sw + c(3) * s2w));
    den = (1 + c(4) * cw + c(5) * c2w) ...
          + 1i * (-(c(4) * sw + c(5) * s2w));
    h = num / den;
end

function g = fce_biquad_peak_gain(c, grid)
    best = 0;
    for i = 0:grid
        w = pi * i / grid;
        v = abs(fce_eval_biquad(c, w));
        if v > best
            best = v;
        end
    end
    g = best;
end

% ======================================================================
% Kaiser
% ======================================================================
function beta = fce_kaiser_beta(att)
    a = abs(att);
    if a > 50
        beta = 0.1102 * (a - 8.7);
    elseif a > 21
        beta = 0.5842 * (a - 21) ^ 0.4 + 0.07886 * (a - 21);
    else
        beta = 0;
    end
end

function taps = fce_kaiser_taps(att, trans, fs)
    a = abs(att);
    dw = 2 * pi * trans / fs;
    n = (a - 7.95) / (2.285 * dw) + 1;
    if n < 1, n = 1; end
    taps = ceil(n);
    if mod(taps, 2) == 0
        taps = taps + 1;
    end
end

% ======================================================================
% FIR (port of fce_fir.c)
% ======================================================================
function w = fce_window_value(win, n, N, beta, alpha)
    if strcmp(win, 'rectangular')
        w = 1;
    elseif strcmp(win, 'hann')
        w = 0.5 * (1 - cos(2 * pi * n / (N - 1)));
    elseif strcmp(win, 'hamming')
        w = 0.54 - 0.46 * cos(2 * pi * n / (N - 1));
    elseif strcmp(win, 'blackman')
        w = 0.42 - 0.5 * cos(2 * pi * n / (N - 1)) ...
            + 0.08 * cos(4 * pi * n / (N - 1));
    elseif strcmp(win, 'kaiser')
        t = 2 * n / (N - 1) - 1;
        w = fce_i0(beta * sqrt(1 - t * t)) / fce_i0(beta);
    elseif strcmp(win, 'blackman_harris')
        w = 0.35875 - 0.48829 * cos(2 * pi * n / (N - 1)) ...
            + 0.14128 * cos(4 * pi * n / (N - 1)) ...
            - 0.01168 * cos(6 * pi * n / (N - 1));
    elseif strcmp(win, 'bartlett')
        x = 2 * n / (N - 1) - 1;
        w = 1 - abs(x);
    elseif strcmp(win, 'tukey')
        nn = n; nmax = N - 1;
        if alpha <= 0
            w = 1;
        elseif alpha >= 1
            w = 0.5 * (1 - cos(2 * pi * nn / nmax));
        elseif nn < alpha * nmax * 0.5
            w = 0.5 * (1 + cos(pi * (2 * nn / (alpha * nmax) - 1)));
        elseif nn >= nmax * (1 - alpha * 0.5)
            w = 0.5 * (1 + cos(pi * (2 * nn / (alpha * nmax) ...
                                     - 2 / alpha + 1)));
        else
            w = 1;
        end
    else
        error('FilterCoeff:spec', 'unknown window');
    end
end

function h = fce_fir_ideal(type, fs, fc1, fc2, N)
    M = 0.5 * (N - 1);
    h = zeros(1, N);
    for n = 0:N-1
        m = n - M;
        v = 0;
        if strcmp(type, 'lowpass')
            v = (2 * fc1 / fs) * fce_sinc(2 * fc1 * m / fs);
        elseif strcmp(type, 'highpass')
            v = -(2 * fc1 / fs) * fce_sinc(2 * fc1 * m / fs);
            if m == 0, v = v + 1; end
        elseif strcmp(type, 'bandpass')
            v = (2 * fc2 / fs) * fce_sinc(2 * fc2 * m / fs) ...
                - (2 * fc1 / fs) * fce_sinc(2 * fc1 * m / fs);
        elseif strcmp(type, 'bandstop')
            v = (2 * fc1 / fs) * fce_sinc(2 * fc1 * m / fs) ...
                - (2 * fc2 / fs) * fce_sinc(2 * fc2 * m / fs);
            if m == 0, v = v + 1; end
        elseif strcmp(type, 'hilbert')
            if fc1 > 0 && fc2 > fc1
                w1 = 2 * pi * fc1 / fs; w2 = 2 * pi * fc2 / fs;
            else
                w1 = 0; w2 = pi;
            end
            if m == 0
                v = 0;
            else
                v = (cos(w1 * m) - cos(w2 * m)) / (pi * m);
            end
        elseif strcmp(type, 'differentiator')
            if fc1 > 0 && fc2 > fc1
                w1 = 2 * pi * fc1 / fs; w2 = 2 * pi * fc2 / fs;
            else
                w1 = 0; w2 = pi;
            end
            if m == 0
                v = 0;
            else
                v = (w2 * cos(w2 * m) - w1 * cos(w1 * m) ...
                     - (sin(w2 * m) - sin(w1 * m)) / m) / (pi * m);
            end
        end
        h(n + 1) = v;
    end
end

function g = fce_fir_gain_at(h, w)
    re = 0; im = 0;
    for n = 1:numel(h)
        a = w * (n - 1);
        re = re + h(n) * cos(a);
        im = im - h(n) * sin(a);
    end
    g = sqrt(re * re + im * im);
end

function g = fce_fir_peak_in_band(h, fs, f_lo, f_hi)
    grid = 256;
    if f_lo <= 0, f_lo = 0; end
    if f_hi >= 0.5 * fs, f_hi = 0.5 * fs; end
    if f_hi <= f_lo, f_hi = 0.5 * fs; end
    best = 0; f_best = 0;
    for i = 0:grid
        f = f_lo + (f_hi - f_lo) * i / grid;
        gv = fce_fir_gain_at(h, 2 * pi * f / fs);
        if gv > best
            best = gv; f_best = f;
        end
    end
    gr = 0.6180339887498948;
    a = f_best - (f_hi - f_lo) / grid;
    b = f_best + (f_hi - f_lo) / grid;
    if a < f_lo, a = f_lo; end
    if b > f_hi, b = f_hi; end
    c = b - gr * (b - a);
    d = a + gr * (b - a);
    fc = fce_fir_gain_at(h, 2 * pi * c / fs);
    fd = fce_fir_gain_at(h, 2 * pi * d / fs);
    for i = 1:60
        if fc > fd
            b = d; d = c; fd = fc;
            c = b - gr * (b - a);
            fc = fce_fir_gain_at(h, 2 * pi * c / fs);
        else
            a = c; c = d; fc = fd;
            d = a + gr * (b - a);
            fd = fce_fir_gain_at(h, 2 * pi * d / fs);
        end
        if b - a < 1e-12 * (1 + abs(b))
            break;
        end
    end
    g = 0.5 * (fc + fd);
end

function r = fce_fir_design(s)
    fs = s.fs;
    type = s.fir_type;
    window = s.window;
    fc1 = s.fc1;
    fc2 = s.fc2;
    N = s.num_taps;
    beta = 0;
    if strcmp(window, 'kaiser')
        if s.kaiser_beta > 0
            beta = s.kaiser_beta;
        elseif s.stopband_atten_db > 0
            beta = fce_kaiser_beta(s.stopband_atten_db);
        end
    end
    ideal = fce_fir_ideal(type, fs, fc1, fc2, N);
    win = zeros(1, N);
    for n = 1:N
        win(n) = fce_window_value(window, n - 1, N, beta, s.tukey_alpha);
    end
    h = ideal .* win;

    norm = s.normalization;
    if strcmp(norm, 'auto')
        if strcmp(type, 'lowpass') || strcmp(type, 'bandstop')
            norm = 'dc';
        elseif strcmp(type, 'highpass') || strcmp(type, 'differentiator')
            norm = 'nyquist';
        elseif strcmp(type, 'bandpass') || strcmp(type, 'hilbert')
            norm = 'passband_peak';
        else
            norm = 'dc';
        end
    end

    if strcmp(norm, 'dc')
        norm_gain = fce_fir_gain_at(h, 0);
    elseif strcmp(norm, 'nyquist')
        norm_gain = fce_fir_gain_at(h, pi);
    elseif strcmp(norm, 'passband_peak')
        if strcmp(type, 'bandpass') || strcmp(type, 'bandstop')
            f_lo = fc1; f_hi = fc2;
        elseif strcmp(type, 'hilbert')
            if fc1 > 0 && fc2 > fc1
                f_lo = fc1; f_hi = fc2;
            else
                f_lo = 0; f_hi = 0.5 * fs;
            end
        else
            f_lo = 0; f_hi = 0.5 * fs;
        end
        norm_gain = fce_fir_peak_in_band(h, fs, f_lo, f_hi);
    else  % 'none'
        norm_gain = 1;
    end

    if ~(norm_gain > 0) || ~(norm_gain < 1e300)
        error('FilterCoeff:numerical', 'normalization failed');
    end
    h = h / norm_gain;

    r.h = h;
    r.num_taps = N;
    r.kaiser_beta = beta;
    r.norm_factor = 1 / norm_gain;
    r.normalization = norm;
    if strcmp(type, 'hilbert')
        r.symmetry = 'III';
    elseif strcmp(type, 'differentiator')
        r.symmetry = 'IV';
    elseif mod(N, 2) == 1
        r.symmetry = 'I';
    else
        r.symmetry = 'II';
    end
end

% ======================================================================
% IIR prototypes (port of fce_iir.c)
% ======================================================================
function [p, k] = fce_proto_butter(n)
    p = zeros(1, n);
    for i = 1:n
        m = 2 * (i - 1) - n + 1;
        th = pi * m / (2 * n);
        p(i) = -exp(1i * th);   % -cos(th) - i sin(th)
    end
    k = 1;
end

function [p, k] = fce_proto_cheb1(n, rp)
    eps = sqrt(10 ^ (0.1 * rp) - 1);
    mu = asinh(1 / eps) / n;
    prod = 1;
    p = zeros(1, n);
    for i = 1:n
        m = 2 * (i - 1) - n + 1;
        th = pi * m / (2 * n);
        sh = sinh(mu); ch = cosh(mu);
        p(i) = -sh * cos(th) - 1i * ch * sin(th);
        prod = prod * ((sh * cos(th))^2 + (ch * sin(th))^2);
    end
    k = sqrt(prod);
    if mod(n, 2) == 0
        k = k / sqrt(1 + eps * eps);
    end
end

function [p, z, k] = fce_proto_cheb2(n, rs)
    de = 1 / sqrt(10 ^ (0.1 * rs) - 1);
    mu = asinh(1 / de) / n;
    z = [];
    p = zeros(1, n);
    if mod(n, 2) == 1
        nz_total = n - 1;
    else
        nz_total = n;
    end
    for i = 0:nz_total-1
        if mod(n, 2) == 1
            half = (n - 1) / 2;
            if i < half
                m = -(2 * (half - i));
            else
                m = 2 * (i - half + 1);
            end
        else
            m = 2 * i - n + 1;
        end
        th = pi * m / (2 * n);
        z = [z, 1i / sin(th)];
    end
    for i = 1:n
        m = 2 * (i - 1) - n + 1;
        th = pi * m / (2 * n);
        sh = sinh(mu); ch = cosh(mu);
        p(i) = -1 / (sh * cos(th) + 1i * ch * sin(th));
    end
    prod_p = 1; prod_z = 1;
    for i = 1:n
        prod_p = prod_p * (-p(i));
    end
    for i = 1:numel(z)
        prod_z = prod_z * (-z(i));
    end
    k = real(prod_p / prod_z);
end

function [p, z, k] = fce_proto_ellip(n, rp, rs)
    if n == 1
        p0 = -sqrt(1 / (10 ^ (0.1 * rp) - 1));
        p = p0; z = []; k = -p0;
        return;
    end
    eps_sq = 10 ^ (0.1 * rp) - 1;
    eps = sqrt(eps_sq);
    ck1_sq = eps_sq / (10 ^ (0.1 * rs) - 1);
    if ~(ck1_sq > 0) || ck1_sq >= 1
        error('FilterCoeff:spec', 'invalid elliptic spec');
    end
    val0 = fce_ellipk(ck1_sq);
    m = fce_ellipdeg(n, ck1_sq);
    if ~(m > 0) || ~(m < 1)
        error('FilterCoeff:numerical', 'elliptic degree failed');
    end
    capk = fce_ellipk(m);
    z = [];
    for j = 0:((n+1)/2 - 1)
        if mod(n, 2) == 1
            jv = 2 * j;
        else
            jv = 2 * j + 1;
        end
        u = jv * capk / n;
        [sn, ~, ~] = fce_ellipj(u, m);
        if abs(sn) > 1e-14
            z = [z, 1i / (sqrt(m) * sn)];
        end
    end
    cnt = numel(z);
    for j = 1:cnt
        z = [z, conj(z(j))];
    end
    w = fce_arc_jac_sn(1i / eps, ck1_sq);
    if ~(isfinite(real(w)) && isfinite(imag(w))) || abs(real(w)) > 1e-9
        error('FilterCoeff:numerical', 'arc_jac failed');
    end
    r = imag(w);
    v0 = capk * r / (n * val0);
    [sv, cv, dv] = fce_ellipj(v0, 1 - m);
    p = [];
    for j = 0:((n+1)/2 - 1)
        if mod(n, 2) == 1
            jv = 2 * j;
        else
            jv = 2 * j + 1;
        end
        u = jv * capk / n;
        [sn, cn, dn] = fce_ellipj(u, m);
        num = -(cn * dn * sv * cv) - 1i * (sn * dv);
        den = 1 - (dn * sv)^2;
        p = [p, num / den];
    end
    if mod(n, 2) == 1
        norm2 = sum(abs(p).^2);
        keep = [];
        for j = 1:numel(p)
            if abs(imag(p(j))) > 1e-14 * sqrt(norm2)
                keep = [keep, p(j)];
            end
        end
        for j = 1:numel(keep)
            p = [p, conj(keep(j))];
        end
    else
        base = numel(p);
        for j = 1:base
            p = [p, conj(p(j))];
        end
    end
    prod_p = 1; prod_z = 1;
    for j = 1:numel(p)
        prod_p = prod_p * (-p(j));
    end
    for j = 1:numel(z)
        prod_z = prod_z * (-z(j));
    end
    k = real(prod_p / prod_z);
    if mod(n, 2) == 0
        k = k / sqrt(1 + eps_sq);
    end
end

function [p, k] = fce_proto_bessel(n)
    p = fce_bessel_poles(n);
    a_last = 1;
    for i = 0:n-1
        a_last = a_last * (2 * n - i);
    end
    for i = 0:n-1
        a_last = a_last / 2;
    end
    target = 1 / sqrt(2);
    lo = 0.1; hi = 1.5;
    for it = 1:200
        prod = 1;
        for i = 1:n
            prod = prod * (1i * hi - p(i));
        end
        if a_last / abs(prod) < target
            break;
        end
        hi = hi * 2;
    end
    for it = 1:200
        mid = 0.5 * (lo + hi);
        prod_lo = 1;
        for i = 1:n
            prod_lo = prod_lo * (1i * lo - p(i));
        end
        g_lo = a_last / abs(prod_lo);
        prod_mid = 1;
        for i = 1:n
            prod_mid = prod_mid * (1i * mid - p(i));
        end
        g_mid = a_last / abs(prod_mid);
        if (g_lo - target) * (g_mid - target) <= 0
            hi = mid;
        else
            lo = mid;
        end
        if hi - lo < 1e-14 * hi
            break;
        end
    end
    nf = 0.5 * (lo + hi);
    p = p / nf;
    k = a_last * nf ^ (-n);
end

% ======================================================================
% analog frequency transformations
% ======================================================================
function wo = fce_prewarp(f, fs)
    wo = 2 * fs * tan(pi * f / fs);
end

function [z, p, k] = fce_tr_lp2lp(z, p, k, wo)
    nz = numel(z); np = numel(p);
    z = z * wo;
    p = p * wo;
    k = k * wo ^ (np - nz);
end

function [z, p, k] = fce_tr_lp2hp(z, p, k, wo)
    nz = numel(z); np = numel(p);
    deg = np - nz;
    pz = 1; pp = 1;
    for i = 1:nz, pz = pz * (-z(i)); end
    for i = 1:np, pp = pp * (-p(i)); end
    k = k * real(pz / pp);
    z = wo ./ z;
    p = wo ./ p;
    z = [z, zeros(1, deg)];
end

function [zo, po, ko] = fce_tr_lp2bp(z, p, k, wo, bw)
    nz = numel(z); np = numel(p);
    deg = np - nz;
    zo = []; po = [];
    for i = 1:nz
        a = z(i) * 0.5 * bw;
        s = sqrt(a * a - wo * wo);
        zo = [zo, a + s];
    end
    for i = 1:nz
        a = z(i) * 0.5 * bw;
        s = sqrt(a * a - wo * wo);
        zo = [zo, a - s];
    end
    for i = 1:np
        a = p(i) * 0.5 * bw;
        s = sqrt(a * a - wo * wo);
        po = [po, a + s];
    end
    for i = 1:np
        a = p(i) * 0.5 * bw;
        s = sqrt(a * a - wo * wo);
        po = [po, a - s];
    end
    zo = [zo, zeros(1, deg)];
    ko = k * bw ^ deg;
end

function [zo, po, ko] = fce_tr_lp2bs(z, p, k, wo, bw)
    nz = numel(z); np = numel(p);
    deg = np - nz;
    zo = []; po = [];
    for i = 1:nz
        a = (0.5 * bw) / z(i);
        s = sqrt(a * a - wo * wo);
        zo = [zo, a + s];
    end
    for i = 1:nz
        a = (0.5 * bw) / z(i);
        s = sqrt(a * a - wo * wo);
        zo = [zo, a - s];
    end
    for i = 1:np
        a = (0.5 * bw) / p(i);
        s = sqrt(a * a - wo * wo);
        po = [po, a + s];
    end
    for i = 1:np
        a = (0.5 * bw) / p(i);
        s = sqrt(a * a - wo * wo);
        po = [po, a - s];
    end
    pz = 1; pp = 1;
    for i = 1:nz, pz = pz * (-z(i)); end
    for i = 1:np, pp = pp * (-p(i)); end
    ko = k * real(pz / pp);
    for i = 1:deg
        zo = [zo, 1i * wo, -1i * wo];
    end
end

function [z, p, k] = fce_tr_bilinear(z, p, k, fs)
    nz = numel(z); np = numel(p);
    deg = np - nz;
    fs2 = 2 * fs;
    pz = 1; pp = 1;
    for i = 1:nz, pz = pz * (fs2 - z(i)); end
    for i = 1:np, pp = pp * (fs2 - p(i)); end
    k = k * real(pz / pp);
    z = (fs2 + z) ./ (fs2 - z);
    p = (fs2 + p) ./ (fs2 - p);
    z = [z, -ones(1, deg)];
end

% ======================================================================
% zpk -> SOS (nearest pairing, scipy convention)
% ======================================================================
function tf = fce_is_real(c)
    tf = abs(imag(c)) <= 1e-12 * (1 + abs(real(c)));
end

function d = fce_dist2(a, b)
    d = (real(a) - real(b))^2 + (imag(a) - imag(b))^2;
end

function consumed = fce_consume(arr, used, n, v)
    consumed = 0;
    need_conj = ~fce_is_real(v);
    if need_conj
        target = 2;
    else
        target = 1;
    end
    for i = 1:n
        if consumed >= target
            break;
        end
        if used(i)
            continue;
        end
        if fce_dist2(arr(i), v) < 1e-20
            used(i) = 1;
            consumed = consumed + 1;
        elseif need_conj && fce_dist2(arr(i), conj(v)) < 1e-20
            used(i) = 1;
            consumed = consumed + 1;
        end
    end
end

function best = fce_nearest_zero(z, used, nz, p1, kind)
    best = nz;
    best_d = 0;
    for i = 1:nz
        if used(i)
            continue;
        end
        if kind == 1 && ~fce_is_real(z(i)), continue; end
        if kind == 2 && fce_is_real(z(i)), continue; end
        d = fce_dist2(z(i), p1);
        if best == nz || d < best_d
            best = i;
            best_d = d;
        end
    end
end

function sos = fce_section_make(zs, ps)
    b = [0, 0, 1];
    a = [0, 0, 1];
    if numel(zs) == 1
        b(1) = 0; b(2) = 1; b(3) = -real(zs(1));
    elseif numel(zs) == 2
        b(1) = 1;
        b(2) = -(real(zs(1)) + real(zs(2)));
        b(3) = real(zs(1)) * real(zs(2)) - imag(zs(1)) * imag(zs(2));
    end
    if numel(ps) == 1
        a(1) = 0; a(2) = 1; a(3) = -real(ps(1));
    elseif numel(ps) == 2
        a(1) = 1;
        a(2) = -(real(ps(1)) + real(ps(2)));
        a(3) = real(ps(1)) * real(ps(2)) - imag(ps(1)) * imag(ps(2));
    end
    sos = [b(1), b(2), b(3), a(2), a(3)];
end

function [sos, k] = fce_zpk2sos(z, p, k)
    dp = p;
    dz = z;
    npad = numel(dp);
    nzpad = numel(dz);
    while nzpad < npad
        dz = [dz, 0 + 0i];
        nzpad = nzpad + 1;
    end
    while npad < nzpad
        dp = [dp, 0 + 0i];
        npad = npad + 1;
    end
    if mod(npad, 2) == 1
        dp = [dp, 0 + 0i];
        dz = [dz, 0 + 0i];
        npad = npad + 1;
        nzpad = nzpad + 1;
    end
    ns = (npad + 1) / 2;
    up = zeros(1, npad);
    uz = zeros(1, nzpad);
    np_rem = npad;
    nz_rem = nzpad;
    sos = zeros(ns, 5);
    for si = ns:-1:1
        p1i = npad;
        worst = 1e300;
        for j = 1:npad
            if up(j), continue; end
            d = abs(1 - abs(dp(j)));
            if p1i == npad || d < worst
                p1i = j; worst = d;
            end
        end
        if p1i == npad
            error('FilterCoeff:numerical', 'no remaining pole');
        end
        p1 = dp(p1i);
        np_rem = np_rem - fce_consume(dp, up, npad, p1);
        nreal_p = 0;
        for j = 1:npad
            if ~up(j) && fce_is_real(dp(j))
                nreal_p = nreal_p + 1;
            end
        end
        nreal_z = 0;
        for j = 1:nzpad
            if ~uz(j) && fce_is_real(dz(j))
                nreal_z = nreal_z + 1;
            end
        end
        if fce_is_real(p1) && nreal_p == 0
            z1i = fce_nearest_zero(dz, uz, nzpad, p1, 1);
            if z1i < nzpad
                nz_rem = nz_rem - fce_consume(dz, uz, nzpad, dz(z1i));
            end
            if z1i < nzpad
                zs = [dz(z1i), 0 + 0i];
            else
                zs = [0 + 0i, 0 + 0i];
            end
            ps = [p1, 0 + 0i];
            sos(si, :) = fce_section_make(zs, ps);
        elseif (np_rem + 1 == nz_rem) && ~fce_is_real(p1) ...
                && nreal_p == 1 && nreal_z == 1
            z1i = fce_nearest_zero(dz, uz, nzpad, p1, 2);
            if z1i == nzpad
                error('FilterCoeff:numerical', 'pairing failed');
            end
            nz_rem = nz_rem - fce_consume(dz, uz, nzpad, dz(z1i));
            zs = [dz(z1i), conj(dz(z1i))];
            ps = [p1, conj(p1)];
            sos(si, :) = fce_section_make(zs, ps);
        else
            if fce_is_real(p1)
                p2i = npad;
                worst2 = 1e300;
                for j = 1:npad
                    if up(j) || ~fce_is_real(dp(j)), continue; end
                    d = abs(1 - abs(dp(j)));
                    if p2i == npad || d < worst2
                        p2i = j; worst2 = d;
                    end
                end
                if p2i == npad
                    error('FilterCoeff:numerical', 'no real pole 2');
                end
                p2 = dp(p2i);
                np_rem = np_rem - fce_consume(dp, up, npad, p2);
            else
                p2 = conj(p1);
            end
            z1i = fce_nearest_zero(dz, uz, nzpad, p1, 0);
            if z1i == nzpad
                sos(si, :) = fce_section_make([], [p1, p2]);
                continue;
            end
            z1 = dz(z1i);
            nz_rem = nz_rem - fce_consume(dz, uz, nzpad, z1);
            if ~fce_is_real(z1)
                zs = [z1, conj(z1)];
            else
                zs = [z1];
                z2i = fce_nearest_zero(dz, uz, nzpad, p1, 1);
                if z2i < nzpad
                    zs = [zs, dz(z2i)];
                    nz_rem = nz_rem - fce_consume(dz, uz, nzpad, dz(z2i));
                end
            end
            sos(si, :) = fce_section_make(zs, [p1, p2]);
        end
    end
end

% ======================================================================
% IIR auto order + design
% ======================================================================
function a = fce_iir_auto_order(s)
    fs = s.fs;
    type = s.iir_type;
    fam = s.iir_family;
    gpass_db = s.passband_ripple_db;
    if ~(gpass_db > 0), gpass_db = 3; end
    gstop_db = s.stopband_atten_db;
    if ~(gstop_db > gpass_db)
        error('FilterCoeff:spec', 'need gstop > gpass');
    end

    if strcmp(type, 'lowpass')
        wp1 = s.fc1; ws1 = s.edge1_hz;
        if ~(ws1 > wp1), error('FilterCoeff:spec', 'edges'); end
        passb1 = fce_prewarp(wp1, fs);
        stopb1 = fce_prewarp(ws1, fs);
        nat = stopb1 / passb1;
        df1 = 0; df2 = 0;
    elseif strcmp(type, 'highpass')
        wp1 = s.fc1; ws1 = s.edge1_hz;
        if ~(ws1 < wp1), error('FilterCoeff:spec', 'edges'); end
        passb1 = fce_prewarp(wp1, fs);
        stopb1 = fce_prewarp(ws1, fs);
        nat = passb1 / stopb1;
        df1 = 0; df2 = 0;
    elseif strcmp(type, 'bandpass')
        wp1 = s.fc1; wp2 = s.fc2;
        ws1 = s.edge1_hz; ws2 = s.edge2_hz;
        if ~(ws1 < wp1 && wp2 < ws2), error('FilterCoeff:spec', 'edges'); end
        passb1 = fce_prewarp(wp1, fs);
        passb2 = fce_prewarp(wp2, fs);
        stopb1 = fce_prewarp(ws1, fs);
        stopb2 = fce_prewarp(ws2, fs);
        n1 = abs((stopb1 * stopb1 - passb1 * passb2) / ...
                 (stopb1 * (passb1 - passb2)));
        n2 = abs((stopb2 * stopb2 - passb1 * passb2) / ...
                 (stopb2 * (passb1 - passb2)));
        nat = min(n1, n2);
        df1 = 0; df2 = 0;
    elseif strcmp(type, 'bandstop')
        ws1 = s.fc1; ws2 = s.fc2;
        wp1 = s.edge1_hz; wp2 = s.edge2_hz;
        if ~(wp1 < ws1 && ws2 < wp2), error('FilterCoeff:spec', 'edges'); end
        passb1 = fce_prewarp(wp1, fs);
        passb2 = fce_prewarp(wp2, fs);
        stopb1 = fce_prewarp(ws1, fs);
        stopb2 = fce_prewarp(ws2, fs);
        n1 = abs((stopb1 * (passb1 - passb2)) / ...
                 (stopb1 * stopb1 - passb1 * passb2));
        n2 = abs((stopb2 * (passb1 - passb2)) / ...
                 (stopb2 * stopb2 - passb1 * passb2));
        nat = min(n1, n2);
        df1 = 0; df2 = 0;
    else
        error('FilterCoeff:spec', 'bad iir_type');
    end
    if ~(nat > 1), error('FilterCoeff:spec', 'bad nat'); end
    gpass = 10 ^ (0.1 * gpass_db);
    gstop = 10 ^ (0.1 * gstop_db);
    if strcmp(fam, 'butterworth')
        n = ceil(log10((gstop - 1) / (gpass - 1)) / (2 * log10(nat)));
    elseif strcmp(fam, 'chebyshev1') || strcmp(fam, 'chebyshev2')
        n = ceil(fce_acosh(sqrt((gstop - 1) / (gpass - 1))) / fce_acosh(nat));
    elseif strcmp(fam, 'elliptic')
        arg1_sq = (gpass - 1) / (gstop - 1);
        arg0 = 1 / nat;
        d00 = fce_ellipk(arg0 * arg0);
        d01 = fce_ellipkm1(arg0 * arg0);
        d10 = fce_ellipk(arg1_sq);
        d11 = fce_ellipkm1(arg1_sq);
        n = ceil(d00 * d11 / (d01 * d10));
    elseif strcmp(fam, 'bessel')
        error('FilterCoeff:spec', 'bessel has no auto-order');
    else
        error('FilterCoeff:spec', 'bad iir_family');
    end
    if n < 1, error('FilterCoeff:spec', 'bad order'); end
    if n > 24, n = 24; end   % FCE_MAX_AUTO_ORDER

    w0 = 1;
    if strcmp(fam, 'butterworth')
        w0 = (gpass - 1) ^ (-1 / (2 * n));
    end

    if strcmp(type, 'lowpass')
        if strcmp(fam, 'chebyshev2')
            df1 = s.edge1_hz;
        else
            wc = passb1 * w0;
            df1 = fs * atan(wc / (2 * fs)) / pi;
        end
    elseif strcmp(type, 'highpass')
        if strcmp(fam, 'chebyshev2')
            df1 = s.edge1_hz;
        else
            wc = passb1 / w0;
            df1 = fs * atan(wc / (2 * fs)) / pi;
        end
    elseif strcmp(type, 'bandpass')
        bw = passb2 - passb1;
        if strcmp(fam, 'chebyshev2'), w0b = 1; else, w0b = w0; end
        wn_hi = sqrt(w0b * w0b * bw * bw / 4 + passb1 * passb2) + w0b * bw * 0.5;
        wn_lo = sqrt(w0b * w0b * bw * bw / 4 + passb1 * passb2) - w0b * bw * 0.5;
        if strcmp(fam, 'chebyshev2')
            df1 = s.edge1_hz; df2 = s.edge2_hz;
        else
            df1 = fs * atan(wn_lo / (2 * fs)) / pi;
            df2 = fs * atan(wn_hi / (2 * fs)) / pi;
        end
    elseif strcmp(type, 'bandstop')
        bw = passb2 - passb1;
        if strcmp(fam, 'chebyshev2'), w0b = 1; else, w0b = w0; end
        discr = sqrt(bw * bw + 4 * w0b * w0b * passb1 * passb2);
        wn_hi = (bw + discr) / (2 * w0b);
        wn_lo = (discr - bw) / (2 * w0b);
        if strcmp(fam, 'chebyshev2')
            df1 = s.fc1; df2 = s.fc2;
        else
            df1 = fs * atan(wn_lo / (2 * fs)) / pi;
            df2 = fs * atan(wn_hi / (2 * fs)) / pi;
        end
    end
    a.order = n;
    a.design_fc1 = df1;
    a.design_fc2 = df2;
end

function y = fce_acosh(x)
    y = log(x + sqrt(x * x - 1));
end

function r = fce_iir_design(s)
    fs = s.fs;
    fam = s.iir_family;
    type = s.iir_type;
    n = s.order;

    if strcmp(fam, 'butterworth')
        [pp, k] = fce_proto_butter(n);
        pz = [];
    elseif strcmp(fam, 'chebyshev1')
        [pp, k] = fce_proto_cheb1(n, s.passband_ripple_db);
        pz = [];
    elseif strcmp(fam, 'chebyshev2')
        [pp, pz, k] = fce_proto_cheb2(n, s.stopband_atten_db);
    elseif strcmp(fam, 'elliptic')
        [pp, pz, k] = fce_proto_ellip(n, s.passband_ripple_db, ...
                                      s.stopband_atten_db);
    elseif strcmp(fam, 'bessel')
        [pp, k] = fce_proto_bessel(n);
        pz = [];
    else
        error('FilterCoeff:spec', 'bad iir_family');
    end

    df1 = s.design_fc1;
    df2 = s.design_fc2;
    if strcmp(type, 'lowpass')
        [z, p, k] = fce_tr_lp2lp(pz, pp, k, fce_prewarp(df1, fs));
    elseif strcmp(type, 'highpass')
        [z, p, k] = fce_tr_lp2hp(pz, pp, k, fce_prewarp(df1, fs));
    elseif strcmp(type, 'bandpass')
        w1 = fce_prewarp(df1, fs);
        w2 = fce_prewarp(df2, fs);
        [z, p, k] = fce_tr_lp2bp(pz, pp, k, sqrt(w1 * w2), w2 - w1);
    elseif strcmp(type, 'bandstop')
        w1 = fce_prewarp(df1, fs);
        w2 = fce_prewarp(df2, fs);
        [z, p, k] = fce_tr_lp2bs(pz, pp, k, sqrt(w1 * w2), w2 - w1);
    else
        error('FilterCoeff:spec', 'bad iir_type');
    end
    [z, p, k] = fce_tr_bilinear(z, p, k, fs);
    [sos, k] = fce_zpk2sos(z, p, k);
    ns = size(sos, 1);

    sos_order = s.sos_order;
    peaks = zeros(1, ns);
    for i = 1:ns
        peaks(i) = fce_biquad_peak_gain(sos(i, :), 256);
    end
    if strcmp(sos_order, 'pole_radius_desc') || strcmp(sos_order, 'internal_gain')
        key = zeros(1, ns);
        for i = 1:ns
            if strcmp(sos_order, 'internal_gain')
                key(i) = peaks(i);
            else
                [p1, p2] = fce_biquad_poles(sos(i, 4), sos(i, 5));
                key(i) = max(abs(p1), abs(p2));
            end
        end
        for i = 1:ns
            best = i;
            for j = (i + 1):ns
                if strcmp(sos_order, 'pole_radius_desc')
                    take = key(j) > key(best);
                else
                    take = key(j) < key(best);
                end
                if take, best = j; end
            end
            if best ~= i
                tmp = sos(i, :); sos(i, :) = sos(best, :); sos(best, :) = tmp;
                tk = key(i); key(i) = key(best); key(best) = tk;
                tp = peaks(i); peaks(i) = peaks(best); peaks(best) = tp;
            end
        end
    end
    for i = 1:ns
        if i == 1, g = k; else, g = 1; end
        sos(i, 1:3) = sos(i, 1:3) * g;
        peaks(i) = peaks(i) * g;
    end
    r.sos = sos;
    r.order = n;
    r.num_sections = ns;
    r.design_fc1 = df1;
    r.design_fc2 = df2;
    r.section_gains = peaks;
end

function [p1, p2] = fce_biquad_poles(a1, a2)
    disc = a1 * a1 - 4 * a2;
    if disc >= 0
        s = sqrt(disc);
        p1 = 0.5 * (-a1 + s);
        p2 = 0.5 * (-a1 - s);
    else
        s = sqrt(-disc);
        p1 = -0.5 * a1 + 1i * 0.5 * s;
        p2 = -0.5 * a1 - 1i * 0.5 * s;
    end
end
