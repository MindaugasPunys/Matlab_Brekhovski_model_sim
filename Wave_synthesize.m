function Wave = Wave_synthesize(args, input_signal_param, params, freq)
 
RespSpectrum = Response_Spectrum(args, input_signal_param, params, freq); 

Wave = real(ifft(RespSpectrum));

end
