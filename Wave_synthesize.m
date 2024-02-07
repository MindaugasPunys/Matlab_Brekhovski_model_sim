function Wave=Wave_synthesize(args,input_signal_param,params,freq,time)
 
RespSpectrum=Response_Spectrum(args,input_signal_param,params,freq,time); 

Wave=real(ifft(RespSpectrum));

end
