function filtercoeff_demo
% FILTERCOEFF_DEMO  Concrete input -> output examples for the MATLAB port.
%
%   >> filtercoeff_demo
%
% Runs a few small designs and prints the inputs and outputs so you can see
% exactly how to use filtercoeff.m and filtercoeff_quantize.m.
    clc;

    % ================================================================
    % 1) FIR lowpass, Hamming window, 101 taps, Fs = 48 kHz, Fc = 5 kHz
    % ================================================================
    sp = struct('kind','fir','fir_type','lowpass','fs',48000, ...
                'fc1',5000,'num_taps',101,'window','hamming');
    [h, info] = filtercoeff(sp);
    fprintf('=== FIR lowpass (Hamming) ===\n');
    fprintf('  input : Fs=%.0f Hz, Fc=%.0f Hz, taps=%d\n', ...
            sp.fs, sp.fc1, info.num_taps);
    fprintf('  output: first 8 of %d taps:\n', info.num_taps);
    fprintf('          %s\n', num2str(h(1:8), '%.6f '));
    fprintf('  sum(h) (DC gain before/after norm) ~ %.4f\n\n', sum(h));

    % ================================================================
    % 2) IIR 4th-order Butterworth lowpass, Fs = 48 kHz, Fc = 5 kHz
    % ================================================================
    sp = struct('kind','iir','iir_family','butterworth','iir_type', ...
                'lowpass','fs',48000,'fc1',5000,'order',4);
    [sos, info] = filtercoeff(sp);
    fprintf('=== IIR Butterworth lowpass (order 4) ===\n');
    fprintf('  input : Fs=%.0f Hz, Fc=%.0f Hz, order=%d\n', ...
            sp.fs, sp.fc1, info.order);
    fprintf('  output: %d SOS sections, layout {b0,b1,b2,a1,a2}:\n', ...
            info.num_sections);
    for k = 1:info.num_sections
        fprintf('    %s\n', num2str(sos((k-1)*5 + (1:5)), '%.8f '));
    end
    fprintf('\n');

    % ================================================================
    % 3) IIR 6th-order Elliptic bandpass + fixed-point Q31
    % ================================================================
    sp = struct('kind','iir','iir_family','elliptic','iir_type', ...
                'bandpass','fs',48000,'fc1',2000,'fc2',4000,'order',6, ...
                'passband_ripple_db',1.0,'stopband_atten_db',50);
    [sos, ~] = filtercoeff(sp);
    [q, scale, intbits] = filtercoeff_quantize(sos, 'q31', 'symmetric', ...
                                               numel(sos));
    fprintf('=== IIR Elliptic bandpass (order 6) + Q31 ===\n');
    fprintf('  input : Fs=48k, passband [%.0f,%.0f] Hz, rp=%.1f, rs=%.0f\n', ...
            sp.fc1, sp.fc2, sp.passband_ripple_db, sp.stopband_atten_db);
    fprintf('  output: %d sections (float + Q31)\n', numel(sos)/5);
    fprintf('          Q31 scale=%.4g, int_bits=%d\n', scale, intbits);
    fprintf('          first section Q31: %s\n\n', num2str(q(1:5), '%d '));

    % ================================================================
    % 4) Auto-order (no order given) - Butterworth lowpass
    % ================================================================
    sp = struct('kind','iir','iir_family','butterworth','iir_type', ...
                'lowpass','fs',10000,'fc1',1000,'order',0, ...
                'edge1_hz',2000,'passband_ripple_db',3.0, ...
                'stopband_atten_db',60);
    [sos, info] = filtercoeff(sp);
    fprintf('=== IIR auto-order (Butterworth lowpass) ===\n');
    fprintf('  input : Fs=10k, passband 0..1000 Hz, stopband from 2000 Hz\n');
    fprintf('  output: computed order = %d, %d sections\n', ...
            info.order, info.num_sections);
    fprintf('          (SciPy buttord gives the same order)\n');
end
