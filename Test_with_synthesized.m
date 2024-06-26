% Bandymas aproksimuoti testinius signalus.
clear all
% Settings of simulation (START)
% addpath('F:\_SMART\Modelio_sudarymas_ir_tikrinimas\Matlab_model_funkcijos\'); %Funkciju nuoroda
% addpath('F:\Dropbox\Mokslas\_SMART\FittingAlgs\Scripts\');

SaveToFile=0; % if set to 1 will save data to savefile
DrawFigures1=1; % draw figures if set to 1
savefile='Test_3nV_TempVar02nan.mat';% 'Test_3nV_2048.mat';
loadfile=''; %'Test_3nV_2048.mat';
FixedNoise=0; % if 1 then wave record with noise is loaded from file and remains fixed
global AirArguments;
AirArguments=0;  % if set to 1 Temperature and Pressure are used instead of c1  and ro1

TestCycles_New=20;  % number of times to repeat approximation
% Settings of simulation (END)
%%
% Bendri parametrai
N_points=1024*2;

% https://en.wikipedia.org/wiki/Speed_of_sound
% Temperature=20; % [Celsijus] temperatura
% Pressure=101.325e3; % [Pa] standard atmospheris pressure
% c1=331.38*sqrt(1+Temperature/273.15) %
% gama=1.4; % air adiabatic constant
% Ks=gama*Pressure; % coefficient of stiffness
% ro1=Ks/(c1^2)

F_sampling=10e6;
T_sampling=1/F_sampling;

% Excitation signal
Magnitude=1000;

F1_chirp=0.3e6;  %0.3
F2_chirp=0.95e6;  % 0.95
T_chirp=10e-6;
T_StartShift=20e-6;
T_All=(N_points-1)*T_sampling;
% Receiver electronics input noise
en=0;  %3e-9; % V/sqrt(Hz)

XducF1=F1_chirp;
XducF2=F2_chirp;

[Reference,time]=Excitation_Chirp(F1_chirp,F2_chirp,T_chirp,T_StartShift,Magnitude,XducF1,XducF2,T_sampling,T_All);

% ExcitationPulseWidth=20;
% StartShift=2024;
% Reference=[zeros(1,StartShift), Magnitude*ones(1,ExcitationPulseWidth), -Magnitude*ones(1,ExcitationPulseWidth), zeros(1,N_points-2*ExcitationPulseWidth-StartShift)];
% time=(0:T_sampling:T_sampling*(length(Reference)-1))*1e6; %us

% Synthesized Reference signal
% alfa0=18.382; % atteniuation coefficient limits
% freq0=461426;% normalised (resonant?) frquency limits
% V_sluoksnio=340;% Ultrasound velocity of layer limits
% n=1.9;% exponentiation number of attenuations limits
% ro2=ro1;% Layer density limits
% h=0.001174;% Layer thikness limits
% Original_arguments=[alfa0,freq0,V_sluoksnio,n,ro2,h];
% Wave_Ref=Wave_synthesize(Original_arguments,Excitation, parameters);

%%
% Synthesized Measured signal
alfa0=18.382; % attenuation coefficient
freq0=461426;% normalized (resonant?) frequency
V_sluoksnio=1080;% Ultrasound velocity of layer
n=1.9;% exponentiation number of attenuations
ro2=1124;% Layer density
h=0.001174;% Layer thickness

if (AirArguments==1)
    Temp_n_c=20;
    Press_n_c=101.325e3;
    Press_min=94e3;
    Press_max=105e3;
    Temperature=20; %[-20:10:40]
    Pressure=Press_n_c;
    OriginalArgs=[alfa0,freq0,V_sluoksnio,n,ro2,h,Temperature,Pressure];
    LB_coeff=-[50, 99, 50, 50, 50, 99, 50, 20]; % [procentais]
    UB_coeff=[50, 100, 50, 50, 50, 100, 50, 20]; % [procentais]
    acquisition_parameters=[F_sampling];
    acquisition_parameters_error=[0]; % [procentais]
else
    c1=344;    % First media ultrasound velocity (Air)
    ro1=1.2;   % First media density (Air)  kg/m^3
    OriginalArgs=[alfa0,freq0,V_sluoksnio,n,ro2,h];
    LB_coeff=-[50, 50, 50, 50, 50, 90]; % [procentais]
    UB_coeff=[50, 90, 50, 50, 50, 100]; % [procentais]
    acquisition_parameters=[F_sampling,c1,ro1];
    % Errors of acquisition parameters
    %e_c1=[-20 -10 0 10 20]; % [m/s]
    e_c1=[-10 0 10]; % [m/s]
    e_F_sampling=zeros(size(e_c1)); % sampling frequency is known precisely
    e_ro1=zeros(size(e_c1));
    acquisition_parameters_error=[e_F_sampling/F_sampling*100;e_c1/c1*100;e_ro1/ro1*100]; % [procentais]
end
Num_Parame_Error_Levels=length(e_c1);
%%

LB=OriginalArgs.*(1+LB_coeff/100);
UB=OriginalArgs.*(1+UB_coeff/100);

nfft = 2^nextpow2(length(Reference));
HalfXcorrLength=floor(nfft/2);
nr=(1:HalfXcorrLength+1)-1;
freq_axis(1:HalfXcorrLength+1)=F_sampling/nfft*nr;
nr=(HalfXcorrLength+2:nfft)-(nfft+1);
freq_axis(HalfXcorrLength+2:nfft)=F_sampling/nfft*nr;

Wave_Meas_WithoutNoise=Wave_synthesize(OriginalArgs,Reference,acquisition_parameters,freq_axis);
% Wave_Meas_WithoutNoise2=Wave_synthesizeModel2(OriginalArgs,Reference,acquisition_parameters,freq_axis);
% max(Wave_Meas_WithoutNoise-Wave_Meas_WithoutNoise2)
% figure
% %subplot(211)
% plot(Wave_Meas_WithoutNoise,'r');
% hold on
% %subplot(212)
% plot(Wave_Meas_WithoutNoise2,'b');
% grid
% pause

% GaindB=60;  %dB
% Gain=10^(GaindB/20);
Gain=max(Reference)/max(Wave_Meas_WithoutNoise);

% ======== Settings for fitting function (START) ===================
% search limits descriptions:
% LB(1)=18.382/1.5; UB(1)=18.382*1.5;  % atteniuation coefficient limits
% LB(2)=461426/3; UB(2)=461426*3;      % normalised (resonant?) frquency limits
% LB(3)=1080/1.5; UB(3)=1080*1.5;      % Ultrasound velocity of layer limits
% LB(4)=1.9/1.5; UB(4)=1.9*1.5;        % exponentiation number of attenutions limits
% LB(5)=1124/1.7; UB(5)=1124*1.7;      % Layer density limits
% LB(6)=0.001174/10; UB(6)=0.001174*2; % Layer thickness limits

%LB(1)=-18.382/1.5; UB(1)=18.382*1.5; % Layer thickness limits
nvariables=length(OriginalArgs);

options = optimoptions('particleswarm','SwarmSize',200,'HybridFcn',@fmincon, 'MaxStallIterations',100,'Display','off','FunctionTolerance',1e-6);
%options = optimoptions('particleswarm','SwarmSize',100,'HybridFcn',@fmincon, 'MaxStallIterations',100,'Display','off','FunctionTolerance',1e-8);
%options = optimoptions('particleswarm','SwarmSize',100,'HybridFcn',@patternsearch, 'MaxStallIterations',100,'Display','off','FunctionTolerance',1e-8);

% ======== Settings for fitting function (END) ===================
%%
tstart_cpu = tic;

for m=1:Num_Parame_Error_Levels
    % acq_parameters_shifted with error due to not precissely known values
    acq_parameters_shifted(m,:)=acquisition_parameters.*(1+acquisition_parameters_error(:,m)'/100);

    if (FixedNoise==1)
        load(loadfile);
        Noise_Meas=StrangeRecordsNoise(1).record; % take first strange record for analysis
        Wave_Meas=StrangeRecordsWave(1).record; % take first strange record for analysis

        clear StrangeRecordsNoise StrangeRecordsWave WaveApproximError Args_Errors
    end
    TestCycles=TestCycles_New;
    StrangeRecordsNoise=[];
    StrangeRecordsWave=[];
    StrangeRecordsErrorLevel=-60; %[dB]

    for k=1:TestCycles

        if (FixedNoise~=1)
            Noise_Meas=(en*Gain*sqrt(F_sampling))*randn(size(Wave_Meas_WithoutNoise));
            Wave_Meas=Wave_Meas_WithoutNoise+Noise_Meas;
        end

        % subplot(211)
        % plot(time,Reference,'b');grid;
        % title('Reference');
        % subplot(212)
        % plot(time,Wave_Meas,'r');grid;
        % title('Measured');
        % hold on
        % xlabel('time, us');
        % pause

        % FitBVDmodel = @(Get_parameter)  Model(Get_parameter, Dat_wo_tail, Refs, c1, ro1, Fs); % model description function
        Objective_Function = @(func_args)  Objective_Func(func_args,Reference,Wave_Meas,acq_parameters_shifted(m,:),freq_axis); % model description function

        % grafinis tikslo funkcijos atvaizdavimas
        % figure
        % arg_range=LB:(UB-LB)/100:UB;
        % for i=1:length(arg_range)
        %   Obj_func(i)=Objective_Func(arg_range(i),Reference,Wave_Meas,acquisition_parameters);
        % end
        % plot(arg_range,Obj_func,'.k-');
        %Obj_func_0=Objective_Func(OriginalArgs,Reference,Wave_Meas,acquisition_parameters)
        %pause

        FittedArgs(k,:) = particleswarm(Objective_Function,nvariables,LB,UB,options);

        % Err=Objective_Func(FittedArgs,Reference,acquisition_parameters);
        % disp('Function fitting error');
        % disp(num2str(Err));

        Args_Errors(k,:)=(FittedArgs(k,:)-OriginalArgs)./OriginalArgs*100;
        %disp('Original,  fitted, and error (%)');
        %disp('alfa0   freq0   V_sluoksnio      n   ro2    h  c1   ro1');
        %disp(num2str([OriginalArgs;FittedArgs;Args_Errors],3));

        Wave_Fitted=Wave_synthesize(FittedArgs(k,:),Reference,acquisition_parameters,freq_axis);
        WaveApproximError(k,m)=20*log10(sum((Wave_Fitted-Wave_Meas_WithoutNoise).^2)/sum(Wave_Meas_WithoutNoise.^2));
        WaveApproximErrorNoisy(k,m)=20*log10(sum((Wave_Fitted-Wave_Meas).^2)/sum(Wave_Meas.^2));

        disp(['Test no. ' num2str(k) '. Wave Approximation Error =' num2str(WaveApproximError(k,m)) ' dB']);

        % save strange records for further analysis
        if (WaveApproximError(k,m)>StrangeRecordsErrorLevel)
            noise.error=WaveApproximError(k,m);
            noise.record=Noise_Meas;
            StrangeRecordsNoise=[StrangeRecordsNoise; noise];
            meas.error=WaveApproximError(k,m);
            meas.record=Wave_Meas;
            StrangeRecordsWave=[StrangeRecordsWave; meas];
            disp(['Strange record found. Total strange records ' num2str(length(StrangeRecordsNoise))]);
        end

        % figure
        % p(1)=plot(time,Wave_Meas,'b');grid;
        % hold on
        % p(2)=plot(time,Wave_Fitted,'r');grid;
        % legend(p,['Measured';'Fitted  ']);

    end  %  for k=1:TestCycles

    MedianWAE(m)=median(WaveApproximError(:,m));
    MedianWAENoisy(m)=median(WaveApproximErrorNoisy(:,m));
    Percentage_strange_recs(m)=length(StrangeRecordsNoise)/TestCycles*100;

    MedianArgs_Errors(m,:)=median(Args_Errors);
    MeanArgs_Errors(m,:)=mean(Args_Errors);
    stdArgs_Errors(m,:)=std(Args_Errors);

    disp(['============ m=' num2str(m) ' ============================']);
    disp(['Median of Approximation Error=' num2str(MedianWAE(m)) ' dB']);
    disp(['Median of Approximation Error Noisy=' num2str(MedianWAENoisy(m)) ' dB']);
    disp(['Percentage of strange records ' num2str(Percentage_strange_recs(m)) ' %']);
    disp('======================================================');

    if (SaveToFile==1)
        save(savefile,'WaveApproximError','WaveApproximErrorNoisy','Args_Errors','StrangeRecordsNoise',...
            'StrangeRecordsWave','MedianWAE','Percentage_strange_recs','TestCycles');
    end
    %load('Test_3nV.mat');
    %load('Test_30nV.mat');
    %load('Test_0nV.mat');

    if (DrawFigures1==1)
        figure
        pl(1)=plot(WaveApproximError(:,m),'.k-');
        hold on
        pl(2)=plot(WaveApproximErrorNoisy(:,m),'.r-');
        legend(pl,'Refrence noise free','Reference noisy');
        ylabel('Wave Approx. Error, dB');
        xlabel('Test number')
        grid

        if (AirArguments==1)
            TitleStrings={'alfa0','freq0','c2','n','ro2','h','T','P'};
        else
            TitleStrings={'alfa0','freq0','c2','n','ro2','h'};
        end

        figure
        for j=1:length(TitleStrings)
            subplot(length(TitleStrings)/2,2,j)
            plot(Args_Errors(:,j),'.b-');
            ylabel(['Error, [%]']);
            xlabel('Test num')
            title(['Error of ' TitleStrings{j}]);
            grid
        end % for j=1:length(TitleStrings)
    end % if (DrawFigures1==1)

end  %for m=1:Num_Parame_Error_Levels

%% 
ParNo=2; % c1 parameter number 2
ArgNo=6; % 6 - h parameter
%stdMedian_c1_Errors(m,:)=std(MedianArgs_Errors(:,ArgNo));

figure
p1(1)=plot(acq_parameters_shifted(:,ParNo)',MedianArgs_Errors(:,ArgNo)','or-');
hold on
p1(2)=plot(acq_parameters_shifted(:,ParNo)',MeanArgs_Errors(:,ArgNo)','*m-');
p1(3)=plot(acq_parameters_shifted(:,ParNo)',MedianArgs_Errors(:,ArgNo)'-2*stdArgs_Errors(:,ArgNo)','.r--');
p1(4)=plot(acq_parameters_shifted(:,ParNo)',MedianArgs_Errors(:,ArgNo)'+2*stdArgs_Errors(:,ArgNo)','.r--');
xlabel('c1, m/s');
ylabel('Error h, %');
legend(p1,'Median','Mean','Median-2std','Median+2std');
grid

figure
plot(acq_parameters_shifted(:,ParNo)',MedianArgs_Errors(:,ArgNo)','or-');
xlabel('c1, m/s');
ylabel('Median (Error h), %');
title('Median')
grid
%===============================================================

ArgNo=5; % ro2 - h parameter
%stdMedian_c1_Errors(m,:)=std(MedianArgs_Errors(:,ArgNo));

figure
p1(1)=plot(acq_parameters_shifted(:,ParNo)',MedianArgs_Errors(:,ArgNo)','or-');
hold on
p1(2)=plot(acq_parameters_shifted(:,ParNo)',MeanArgs_Errors(:,ArgNo)','*m-');
p1(3)=plot(acq_parameters_shifted(:,ParNo)',MedianArgs_Errors(:,ArgNo)'-2*stdArgs_Errors(:,ArgNo)','.r--');
p1(4)=plot(acq_parameters_shifted(:,ParNo)',MedianArgs_Errors(:,ArgNo)'+2*stdArgs_Errors(:,ArgNo)','.r--');
xlabel('c1, m/s');
ylabel('Error ro2, %');
legend(p1,'Median','Mean','Median-2std','Median+2std');
grid

figure
plot(acq_parameters_shifted(:,ParNo)',MedianArgs_Errors(:,ArgNo)','or-');
xlabel('c1, m/s');
ylabel('Median (Error ro2), %');
title('Median')
grid

%===============================================================

ArgNo=3; % V_sluoksnio - h parameter
%stdMedian_c1_Errors(m,:)=std(MedianArgs_Errors(:,ArgNo));

figure
p1(1)=plot(acq_parameters_shifted(:,ParNo)',MedianArgs_Errors(:,ArgNo)','or-');
hold on
p1(2)=plot(acq_parameters_shifted(:,ParNo)',MeanArgs_Errors(:,ArgNo)','*m-');
p1(3)=plot(acq_parameters_shifted(:,ParNo)',MedianArgs_Errors(:,ArgNo)'-2*stdArgs_Errors(:,ArgNo)','.r--');
p1(4)=plot(acq_parameters_shifted(:,ParNo)',MedianArgs_Errors(:,ArgNo)'+2*stdArgs_Errors(:,ArgNo)','.r--');
xlabel('c1, m/s');
ylabel('Error c2, %');
legend(p1,'Median','Mean','Median-2std','Median+2std');
grid

figure
plot(acq_parameters_shifted(:,ParNo)',MedianArgs_Errors(:,ArgNo)','or-');
xlabel('c1, m/s');
ylabel('Median (Error c2), %');
title('Median')
grid

telapsed_cpu = toc(tstart_cpu);
disp(['Simulation duration: ' num2str(telapsed_cpu/60) ' minutes']);

