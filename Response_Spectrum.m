function RespSpectrum = Response_Spectrum(args, input_signal_param, params, freq)

nfft = length(freq);
input_fft = fft(input_signal_param, nfft);
h = args(5);
c1 = params(2);
T = TransferFunction(args, params, freq);
RespSpectrum = input_fft .* conj(T) .* exp(1i * 2 * pi * freq * (h / c1));

end
