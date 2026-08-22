function [q, scale, int_bits] = filtercoeff_quantize(coeffs, qformat, strategy, sec_len)
% FILTERCOEFF_QUANTIZE  Fixed-point conversion (port of fce_quant.c).
%
%   [q, scale, int_bits] = filtercoeff_quantize(coeffs, qformat, strategy, sec_len)
%
%   coeffs   : 1xN vector of float coefficients (FIR taps or flat SOS)
%   qformat  : 'q15' or 'q31'
%   strategy : 'symmetric' | 'section_wise' | 'coefficient_wise'
%   sec_len  : coefficients per scaling unit for 'section_wise' (5 for SOS)
%   q        : 1xN integer vector  (q = round(c * scale), saturated)
%   scale    : symmetric scale ([] for section/coefficient-wise)
%   int_bits : effective integer bits
%
% Reconstruction: c_tilde = q / scale  (symmetric) or per-section/coefficient.
    if nargin < 3, strategy = 'symmetric'; end
    if nargin < 4, sec_len = numel(coeffs); end
    n = numel(coeffs);
    if n == 0
        error('FilterCoeff:quant', 'empty coefficient set');
    end
    if strcmp(qformat, 'q15')
        frac = 15;
    elseif strcmp(qformat, 'q31')
        frac = 31;
    else
        error('FilterCoeff:quant', 'qformat must be q15 or q31');
    end
    qmax = 2 ^ frac - 1;
    qmax_d = double(qmax);

    if strcmp(strategy, 'symmetric')
        n_sec = ceil(n / sec_len);
        mx = max(abs(coeffs));
        if mx > 0, scale = qmax_d / mx; else, scale = 1; end
        scales = scale * ones(1, n);
        if mx > 1
            int_bits = ceil(log2(mx));
        else
            int_bits = 0;
        end
    elseif strcmp(strategy, 'section_wise')
        if nargin < 4
            error('FilterCoeff:quant', 'section_wise needs sec_len');
        end
        n_sec = ceil(n / sec_len);
        scales = zeros(1, n);
        scale = [];
        int_bits = 0;
        for sec = 0:n_sec-1
            base = sec * sec_len;
            cnt = min(sec_len, n - base);
            mx = 0;
            for i = 1:cnt
                mx = max(mx, abs(coeffs(base + i)));
            end
            if mx > 0, s = qmax_d / mx; else, s = 1; end
            for i = 1:cnt
                scales(base + i) = s;
            end
        end
    elseif strcmp(strategy, 'coefficient_wise')
        scales = zeros(1, n);
        mx = 0;
        scale = [];
        for i = 1:n
            a = abs(coeffs(i));
            mx = max(mx, a);
            if a > 0
                scales(i) = qmax_d / a;
            else
                scales(i) = 1;
            end
        end
        if mx > 1
            int_bits = ceil(log2(mx));
        else
            int_bits = 0;
        end
    else
        error('FilterCoeff:quant', 'unknown strategy');
    end

    q = zeros(1, n);
    for i = 1:n
        if scales(i) > 0
            qv = fce_round_half(coeffs(i) * scales(i));
        else
            qv = 0;
        end
        qv = max(-qmax, min(qmax, qv));
        q(i) = qv;
    end
end

function y = fce_round_half(x)
    if x >= 0
        y = floor(x + 0.5);
    else
        y = ceil(x - 0.5);
    end
end
