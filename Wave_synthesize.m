function Wave = Wave_synthesize(args, input_signal_param, params, freq)
DEBUG = 0;
 %% Model Transfer function
T = TransferFunction(args, params, freq);

%% Model Response spectrum
c1 = params(2);
h = args(6);

nfft = length(freq);
input_fft = fft(input_signal_param, nfft); % >>> FFT <<<

RespSpectrum = input_fft .* conj(T) .* exp(1i * 2 * pi * freq * (h / c1));

%% Synthesize wave in time domain
Wave = real(ifft(RespSpectrum)); % >>> IFFT <<<

%% TEST
if DEBUG
    % dir = 'C:/GIT/VITIS/Model_sim/solution/csim/report/TransferFunction_csim.log';
    dir = 'C:/GIT/VITIS/Model_sim/solution/csim/report/WaveSynthesis_csim.log';

    % input_fft_hls = plot_xilinx_data(dir, 'input_fft');
    % hold on;
    % plot(real(input_fft), '--')
    % hold on;
    % plot(imag(input_fft), '--')
    % title('input fft'); legend('HLS re', 'HLS im', 'M re', 'M im');
    % xlabel('Indeksas'); ylabel('Amplitudė');
    % 
    % plot_xilinx_data(dir, 'T_conj');
    % hold on;
    % plot(real(conj(T)), '--')
    % hold on;
    % plot(imag(conj(T)), '--')
    % title('T conj'); legend('HLS re', 'HLS im', 'M re', 'M im');
    % xlabel('Indeksas'); ylabel('Amplitudė');

    plot_xilinx_data(dir, 'exp_arg_step_1');
    hold on;
    plot(freq * h, '--');
    title('exp_arg_step_1'); legend('HLS re', 'M re');
    xlabel('Indeksas'); ylabel('Amplitudė');

    plot_xilinx_data(dir, 'exp_arg_step_2');
    hold on;
    plot(freq / c1 * h, '--');
    title('exp_arg_step_2'); legend('HLS re', 'M re');
    xlabel('Indeksas'); ylabel('Amplitudė');

    plot_xilinx_data(dir, 'exp_arg_step_3');
    hold on;
    plot(2 * pi * freq / c1 * h, '--');
    title('exp_arg_step_3'); legend('HLS re', 'M re');
    xlabel('Indeksas'); ylabel('Amplitudė');

    plot_xilinx_data(dir, 'exp_arg');
    hold on;
    plot(2 * pi * freq * (h / c1), '--');
    title('exp arg'); legend('HLS re', 'M re');
    xlabel('Indeksas'); ylabel('Amplitudė');

    plot_xilinx_data(dir, 'exp_sub');
    hold on;
    plot(real(exp(1i * 2 * pi * freq * (h / c1))), '--')
    hold on;
    plot(imag(exp(1i * 2 * pi * freq * (h / c1))), '--')
    title('exp sub'); legend('HLS re', 'HLS im', 'M re', 'M im');
    xlabel('Indeksas'); ylabel('Amplitudė');

    % RespSpectrum_HLS = plot_xilinx_data(dir, 'resp_spectrum');
    % hold on;
    % plot(real(RespSpectrum), '--')
    % hold on;
    % plot(imag(RespSpectrum), '--')
    % title('RespSpectrum'); legend('HLS re', 'HLS im', 'M re', 'M im');
    % xlabel('Indeksas'); ylabel('Amplitudė');
end
end
