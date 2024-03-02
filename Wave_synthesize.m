function Wave = Wave_synthesize(args, input_signal_param, params, freq)
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
end
