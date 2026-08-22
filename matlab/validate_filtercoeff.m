function status = validate_filtercoeff(refdir)
% VALIDATE_FILTERCOEFF  Compare the MATLAB port against the C library.
%
%   validate_filtercoeff            % reference dir defaults to ./reference
%   validate_filtercoeff('/path')   % explicit reference dir
%
% Reads the coefficient CSVs exported from the C library by
% tools/export_matlab_reference.py, re-designs the same filters with the
% MATLAB port (filtercoeff.m) and reports the maximum absolute coefficient
% error for each design. All designs should match to ~1e-12.
%
% status = 0 if every design is within tolerance, 1 otherwise.
    if nargin < 1
        refdir = fullfile(fileparts(mfilename('fullpath')), 'reference');
    end
    tol = 1e-12;

    specs = build_specs();
    ids = fieldnames(specs);
    nfail = 0;
    for i = 1:numel(ids)
        id = ids{i};
        sp = specs.(id);
        csv = fullfile(refdir, [id '.csv']);
        if ~exist(csv, 'file')
            fprintf('SKIP  %-22s (no reference file)\n', id);
            continue;
        end
        ref = dlmread(csv, ',', 1, 0);
        ref = ref(:, 2)';   % coefficient column
        [coeffs, ~] = filtercoeff(sp);
        if numel(coeffs) ~= numel(ref)
            fprintf('FAIL  %-22s length %d vs ref %d\n', id, ...
                    numel(coeffs), numel(ref));
            nfail = nfail + 1;
            continue;
        end
        err = max(abs(coeffs - ref));
        if err < tol
            fprintf('PASS  %-22s max err %.3e\n', id, err);
        else
            fprintf('FAIL  %-22s max err %.3e\n', id, err);
            nfail = nfail + 1;
        end
    end
    fprintf('\nMATLAB vs C reference: %d failures\n', nfail);
    if nfail == 0
        status = 0;
    else
        status = 1;
    end
end

function specs = build_specs()
% Specs mirroring fce_dump.c (tools/reference/fce_dump.c).
    s = struct();

    % ---- FIR ----
    sp = struct('kind','fir','fir_type','lowpass','fs',48000,'fc1',5000, ...
        'num_taps',101,'window','hann','normalization','dc');
    s.fir_lp_hann_odd = sp;

    sp = struct('kind','fir','fir_type','lowpass','fs',48000,'fc1',5000, ...
        'num_taps',101,'window','kaiser','kaiser_beta',7.86, ...
        'normalization','dc');
    s.fir_lp_kaiser = sp;

    sp = struct('kind','fir','fir_type','highpass','fs',48000,'fc1',5000, ...
        'num_taps',101,'window','hamming','normalization','nyquist');
    s.fir_hp_hamming = sp;

    sp = struct('kind','fir','fir_type','bandpass','fs',48000, ...
        'fc1',3000,'fc2',6000,'num_taps',121,'window','hamming', ...
        'normalization','passband_peak');
    s.fir_bp_hamming = sp;

    sp = struct('kind','fir','fir_type','bandstop','fs',48000, ...
        'fc1',3000,'fc2',6000,'num_taps',161,'window','kaiser', ...
        'kaiser_beta',6.2,'normalization','dc');
    s.fir_bs_kaiser = sp;

    sp = struct('kind','fir','fir_type','hilbert','fs',48000, ...
        'fc1',0,'fc2',0,'num_taps',65,'window','hamming', ...
        'normalization','passband_peak');
    s.fir_hilbert = sp;

    sp = struct('kind','fir','fir_type','lowpass','fs',48000,'fc1',5000, ...
        'num_taps',0,'window','kaiser','stopband_atten_db',80, ...
        'transition_hz',1000,'normalization','dc');
    s.fir_lp_kaiser_auto = sp;

    % ---- IIR ----
    sp = struct('kind','iir','iir_family','butterworth','iir_type', ...
        'lowpass','fs',48000,'fc1',5000,'order',4);
    s.iir_butter_lp4 = sp;

    sp = struct('kind','iir','iir_family','butterworth','iir_type', ...
        'highpass','fs',44100,'fc1',5000,'order',8);
    s.iir_butter_hp8 = sp;

    sp = struct('kind','iir','iir_family','butterworth','iir_type', ...
        'bandpass','fs',48000,'fc1',2000,'fc2',4000,'order',4);
    s.iir_butter_bp4 = sp;

    sp = struct('kind','iir','iir_family','chebyshev1','iir_type', ...
        'lowpass','fs',48000,'fc1',5000,'order',5, ...
        'passband_ripple_db',1.0);
    s.iir_cheby1_lp5 = sp;

    sp = struct('kind','iir','iir_family','chebyshev2','iir_type', ...
        'lowpass','fs',48000,'fc1',5000,'order',4, ...
        'stopband_atten_db',60);
    s.iir_cheby2_lp4 = sp;

    sp = struct('kind','iir','iir_family','elliptic','iir_type', ...
        'lowpass','fs',48000,'fc1',5000,'order',4, ...
        'passband_ripple_db',0.5,'stopband_atten_db',60);
    s.iir_ellip_lp4 = sp;

    sp = struct('kind','iir','iir_family','elliptic','iir_type', ...
        'bandpass','fs',48000,'fc1',2000,'fc2',4000,'order',6, ...
        'passband_ripple_db',1.0,'stopband_atten_db',50);
    s.iir_ellip_bp6 = sp;

    sp = struct('kind','iir','iir_family','bessel','iir_type', ...
        'lowpass','fs',48000,'fc1',5000,'order',3);
    s.iir_bessel_lp3 = sp;

    specs = s;
end
