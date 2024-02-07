Objective_Function = @(func_args)  Objective_Func(func_args,Wave_Ref,Wave_Meas, acq_parameters_shifted,freq_axis,time_axis,H_electronics,np1,np2,TRANSFER_MODEL); % model description function
[FittedArgs,approximation_error(i),exitflagPSO,outputPSO] = particleswarm(Objective_Function,nvariables,LB,UB,options);


function error = Objective_Func(args,ref_signal,meas_signal,params,f_axis,time,H_electronics,start_point,end_point,model)
model_signal=Wave_synthesize(args,ref_signal,params,f_axis,time,H_electronics,model);
error = sum( (model_signal(start_point:end_point) - meas_signal(start_point:end_point)).^2)/sum(meas_signal(start_point:end_point).^2);
return

end

function Wave=Wave_synthesize(args,input_signal_param,params,freq,time,H_electronics,model)

RespSpectrum = Response_Spectrum(args, input_signal_param, params, freq, time, H_electronics, model);
Wave=real(ifft(RespSpectrum));

end

function RespSpectrum=Response_Spectrum(args,input_signal_param,params,freq,time,H_electronics,model)
global AirArguments
% spektro skaiciavimai
% nfft = 2^nextpow2(length(input_signal));

nfft=length(freq);
input_fft=fft(input_signal_param,nfft);


h=args(5);
c1=params(2);
T=TransferFunction(args,params,freq,model);
RespSpectrum=input_fft.*conj(T).*exp(1i*2*pi*freq*(h/c1));  %

end

function T=TransferFunction(args,params,freq,model)
global AirArguments

alfa0=args(1);
c2=args(2);
n=args(3);
ro2=args(4);
h=args(5);

fsampl=params(1);
ro1=params(3);
freq0=params(4);
c1=params(2);

alf = alfa0*(abs(freq)./freq0).^n;  % slopinimo koeficientas
C = c2;
k=2*pi*(freq)./C+1i*alf; % randamas kompleksinis banginis skaicius vertinant slopinima
veloc=2*pi*(freq)./real(k);

Z1=c1*ro1; % pirmojo sluoksnio impedansas
Z2=veloc*ro2; % antro sluoksnio impedansas

T=(4*Z1*Z2) ./ ( (Z1+Z2).^2 .*exp(-1i*k*h) - (Z1-Z2).^2.*exp(1i*k*h) ); %Brekhovskikh & Godin 1990 (2.4.18)
T(1)=0;
end