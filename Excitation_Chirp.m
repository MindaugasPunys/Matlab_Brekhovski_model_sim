function [y,t]=Excitation_Chirp(F1,F2,T_width,T_StartShift,Magnitude,XducF1,XducF2,T_Sampling,T_All)
% input data
% F1 chirp start frequency in [Hz]
% F2 chirp end frequency in [Hz]
% T_width - duration of chirp pulse in sec
% T_StartShift - chirp pulse start time in [sec]
% Magnitude - magnitude of chirp pulse
% XducF1 low frequency of transducer
% XducF2 high frequency of transducer
% T_sampling - sampling period in [sec]
% T_All - full signal duration

N_chirp_points=round(T_width/T_Sampling);
t_chirp=0:T_Sampling:T_width;
y_chirp=Magnitude*chirp(t_chirp,F1,T_width,F2,'linear');


N_Points=round(T_All/T_Sampling);
N_start_points=round(T_StartShift/T_Sampling);
t=0:T_Sampling:T_All;
y1=[zeros(1,N_start_points), y_chirp, zeros(1,N_Points-N_start_points-N_chirp_points)];

F_sampling=1/T_Sampling;
% XducParams.f0=4.0e6;
% XducParams.BW_perc=20;
% XducParams.BW=XducParams.f0*XducParams.BW_perc/100;
% XducParams.f1=F1_chirp;  %XducParams.f0-XducParams.BW/2;
% XducParams.f2=F2_chirp;  %XducParams.f0+XducParams.BW/2;
WL=XducF1/(F_sampling/2);
WH=XducF2/(F_sampling/2);
FilterOrder=2;
[Bxduc,Axduc] = butter(FilterOrder,[WL,WH]);
y=filter(Bxduc, Axduc, y1);
end
