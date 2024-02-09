%% SETUP
clear; close all;

%% CONFIGURATION

SaveToFile = 0;

DrawFigures = 1;
figure_num = 1;

savefile = ''; % ''
loadfile = ''; % ''

FixedNoise = 0; % load wave record from loadfile

% 0 - use c1 and ro1
% 1 - use Temperature and Pressure
global AirArguments;
AirArguments = 0;
TestCycles_New = 20;  % times to repeat aproximation

%% CONSTANT VARIABLE DEFINITION

N_points = 1024 * 2;

% https://en.wikipedia.org/wiki/Speed_of_sound
% Temperature = 20; % [degC]
% Pressure = 101.325e3; % [Pa] standard atmospheris pressure
% c1 = 331.38 * sqrt(1 + Temperature / 273.15); % speed of sound
% gama = 1.4; % Air adiabatic constant
% Ks = gama * Pressure; % Coefficient of stiffness
% ro1 = Ks / (c1^2);

F_sampling = 10e6;
T_sampling = 1 / F_sampling;
Magnitude = 1000;

%% EXCITATION CHIRP SIGNAL GENERATION

F1_chirp = 0.3e6;  % 0.3 MHz
F2_chirp = 0.95e6;  % 0.95 MHz
T_chirp = 10e-6; % 10 us
T_StartShift = 20e-6; % 20 us
T_All = (N_points - 1) * T_sampling;
% Receiver electronics input noise
en=0;  %3e-9; % V/sqrt(Hz)

XducF1 = F1_chirp;
XducF2 = F2_chirp;

[Reference, time] = Excitation_Chirp(F1_chirp, F2_chirp, T_chirp, T_StartShift, ...
    Magnitude, XducF1, XducF2, T_sampling, T_All);

if (DrawFigures == 1)
    %% Plot the generated signal
    figure(figure_num); figure_num = figure_num + 1;

    plot(time, Reference);
    xlabel('Time (s)');
    ylabel('Amplitude');
    title('Chirp Excitation Signal');
end


%% SIGNAL PARAMETERS

% Synthesized Measured signal
alfa0 = 18.382; % attenuation coefficient
freq0 = 461426; % normalized (resonant) frequency
V_sluoksnio = 1080; % Ultrasound velocity of layer
n = 1.9; % Exponentiation number of attenuations
ro2 = 1124; % Layer density
h = 0.001174; % Layer thickness

if (AirArguments == 1)
    %% AirArguments == 1; use Temperature and Pressure model
    Temp_n_c = 20;
    Press_n_c = 101.325e3;
    Press_min = 94e3;
    Press_max = 105e3;
    Temperature = 20; %[-20:10:40]
    Pressure = Press_n_c;

    c1 = 331.38 * sqrt(1 + Temperature / 273.15); % speed of sound
    gama = 1.4; % Air adiabatic constant
    Ks = gama * Pressure; % Coefficient of stiffness
    ro1 = Ks / (c1^2);

    OriginalArgs = [alfa0, freq0, V_sluoksnio, n, ro2, h, Temperature, Pressure];
    LB_coeff = -[50, 99, 50, 50, 50, 99, 50, 20]; % [%]
    UB_coeff = [50, 100, 50, 50, 50, 100, 50, 20]; % [%]
    acquisition_parameters = [F_sampling, c1, ro1];
    acquisition_parameters_error = 0; % [%]
    e_c1 = 0;
else
    %% AirArguments == 0; use c1 and ro1 model

    c1 = 344; % First media ultrasound velocity (Air)
    ro1 = 1.2; % First media density (Air)  kg/m^3

    OriginalArgs = [alfa0, freq0, V_sluoksnio, n, ro2, h];
    LB_coeff = -[50, 50, 50, 50, 50, 90]; % [%]
    UB_coeff = [50, 90, 50, 50, 50, 100]; % [%]
    acquisition_parameters = [F_sampling, c1, ro1];
    e_c1 = [-10 0 10]; % [m/s] Errors of acquisition parameters [-20 -10 0 10 20]
    e_F_sampling = zeros(size(e_c1)); % sampling frequency is known precisely
    e_ro1 = zeros(size(e_c1));
    acquisition_parameters_error = [
        e_F_sampling / F_sampling * 100;
        e_c1 / c1 * 100;
        e_ro1 / ro1 * 100]; % [%]
end
Num_Parame_Error_Levels = length(e_c1);

%% MODEL PARAMETERS

LB = OriginalArgs .* (1 + LB_coeff / 100);
UB = OriginalArgs .* (1 + UB_coeff / 100);

nfft = 2 ^ nextpow2(length(Reference));
HalfXcorrLength = floor(nfft / 2);
nr = (1 : HalfXcorrLength+1) - 1;
freq_axis(1 : HalfXcorrLength+1) = F_sampling / nfft * nr;
nr = (HalfXcorrLength+2 : nfft) - (nfft + 1);
freq_axis(HalfXcorrLength+2 : nfft) = F_sampling / nfft * nr;

Wave_Meas_WithoutNoise = Wave_synthesize(OriginalArgs, Reference, acquisition_parameters, freq_axis);


if (DrawFigures == 1)
    %% Plot the generated ultrasound wave
    figure(figure_num); figure_num = figure_num + 1;

    plot(Wave_Meas_WithoutNoise,'r');
    grid
end


% GaindB = 60;  % [dB]
% Gain = 10^(GaindB/20);
Gain = max(Reference) / max(Wave_Meas_WithoutNoise);

%% FITTING FUNCTION SETUP

% Search limits descriptions:
% LB(1) = 18.382 / 1.5;   UB(1) = 18.382 * 1.5;   % atteniuation coefficient limits
% LB(2) = 461426 / 3;     UB(2) = 461426 * 3;     % normalised (resonant) frquency limits
% LB(3) = 1080 / 1.5;     UB(3) = 1080 * 1.5;     % Ultrasound velocity of layer limits
% LB(4) = 1.9 / 1.5;      UB(4) = 1.9 * 1.5;      % exponentiation number of attenutions limits
% LB(5) = 1124 / 1.7;     UB(5) = 1124 * 1.7;     % Layer density limits
% LB(6) = 0.001174 / 10;  UB(6) = 0.001174 * 2;   % Layer thickness limits

%LB(1) = -18.382 / 1.5; UB(1) = 18.382 * 1.5;   % Layer thickness limits
nvariables = length(OriginalArgs);

options = optimoptions('particleswarm', 'SwarmSize', 200, ...
    'HybridFcn', @fmincon, 'MaxStallIterations', 100, ...
    'Display', 'off', 'FunctionTolerance', 1e-6);

% options = optimoptions('particleswarm', 'SwarmSize', 100, ...
%     'HybridFcn', @fmincon, 'MaxStallIterations', 100, ...
%     'Display', 'off', 'FunctionTolerance', 1e-8);

% options = optimoptions('particleswarm', 'SwarmSize', 100, ...
%     'HybridFcn', @patternsearch, 'MaxStallIterations', 100, ...
%     'Display', 'off', 'FunctionTolerance', 1e-8);

%% FITTING PROCESS

tstart_cpu = tic;

for m=1 : Num_Parame_Error_Levels
    %% Model error levels
    % acq_parameters_shifted with error due to not precissely known values
    acq_parameters_shifted(m,:) = acquisition_parameters .* ...
        (1 + acquisition_parameters_error(:,m)' / 100);

    if (FixedNoise == 1)
        load(loadfile);
        Noise_Meas = StrangeRecordsNoise(1).record; % take first strange record for analysis
        Wave_Meas = StrangeRecordsWave(1).record; % take first strange record for analysis
        clear StrangeRecordsNoise StrangeRecordsWave WaveApproximError Args_Errors
    end

    TestCycles = TestCycles_New;
    StrangeRecordsNoise = [];
    StrangeRecordsWave = [];
    StrangeRecordsErrorLevel = -60; %[dB]

    for k = 1:TestCycles
        %% PSO fitting

        if (FixedNoise ~= 1)
            Noise_Meas = (en * Gain * sqrt(F_sampling)) * ...
                randn(size(Wave_Meas_WithoutNoise));
            Wave_Meas = Wave_Meas_WithoutNoise + Noise_Meas;
        end


        if (DrawFigures == 1)
            %% Plot modeled and reference signal
            figure(figure_num); figure_num = figure_num + 1;

            subplot(211)
            plot(time,Reference,'b'); grid on;
            title('Reference');
            subplot(212)
            plot(time,Wave_Meas,'r'); grid on;
            title('Measured');
            hold on
            xlabel('time, us');
            pause(5);
        end

        % FitBVDmodel = @(Get_parameter)  Model(Get_parameter, Dat_wo_tail, Refs, c1, ro1, Fs); % model description function
        Objective_Function = @(func_args) Objective_Func(func_args, Reference, ...
            Wave_Meas, acq_parameters_shifted(m,:), freq_axis); % model description function
        

        if (DrawFigures == 1)
            %% Plot fitness function
            figure(figure_num); figure_num = figure_num + 1;

            arg_range = LB : (UB-LB)/100 : UB;
            for i = 1 : length(arg_range)
                Obj_func(i) = Objective_Func(arg_range(i), Reference, Wave_Meas, acquisition_parameters);
            end
            plot(arg_range, Obj_func, '.k-');
            Obj_func_0 = Objective_Func(OriginalArgs, Reference, Wave_Meas, acquisition_parameters)
            pause
        end

        FittedArgs(k,:) = particleswarm(Objective_Function, nvariables, LB, UB, options); %#ok<*SAGROW>

        % Err=Objective_Func(FittedArgs,Reference,acquisition_parameters);
        % disp('Function fitting error');
        % disp(num2str(Err));

        Args_Errors(k,:) = (FittedArgs(k,:) - OriginalArgs) ./ OriginalArgs * 100;
        %disp('Original,  fitted, and error (%)');
        %disp('alfa0   freq0   V_sluoksnio      n   ro2    h  c1   ro1');
        %disp(num2str([OriginalArgs;FittedArgs;Args_Errors],3));

        Wave_Fitted = Wave_synthesize(FittedArgs(k,:), Reference, acquisition_parameters, freq_axis);
        WaveApproximError(k,m) = 20 * log10(sum((Wave_Fitted - Wave_Meas_WithoutNoise) .^ 2) /...
                                            sum(Wave_Meas_WithoutNoise .^ 2)                );
        WaveApproximErrorNoisy(k,m) = 20 * log10(sum((Wave_Fitted - Wave_Meas) .^ 2) / ...
                                                 sum(Wave_Meas .^ 2));

        disp(['Test no. ' num2str(k) '. Wave Approximation Error =' num2str(WaveApproximError(k,m)) ' dB']);

        % save strange records for further analysis
        if (WaveApproximError(k,m) > StrangeRecordsErrorLevel)
            noise.error = WaveApproximError(k,m);
            noise.record = Noise_Meas;
            StrangeRecordsNoise = [StrangeRecordsNoise; noise]; %#ok<*AGROW>
            meas.error = WaveApproximError(k,m);
            meas.record = Wave_Meas;
            StrangeRecordsWave = [StrangeRecordsWave; meas];
            disp(['Strange record found. Total strange records ' num2str(length(StrangeRecordsNoise))]);
        end

        % figure
        % p(1)=plot(time,Wave_Meas,'b');grid;
        % hold on
        % p(2)=plot(time,Wave_Fitted,'r');grid;
        % legend(p,['Measured';'Fitted  ']);

    end  %  for k=1:TestCycles

    MedianWAE(m) = median(WaveApproximError(:,m));
    MedianWAENoisy(m) = median(WaveApproximErrorNoisy(:,m));
    Percentage_strange_recs(m) = length(StrangeRecordsNoise) / TestCycles * 100;

    MedianArgs_Errors(m,:) = median(Args_Errors);
    MeanArgs_Errors(m,:) = mean(Args_Errors);
    stdArgs_Errors(m,:) = std(Args_Errors);

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

    if (DrawFigures == 1)
        %% Plot measured errors
        figure(figure_num); figure_num = figure_num + 1;
        pl(1) = plot(WaveApproximError(:,m), '.k-');
        hold on
        pl(2) = plot(WaveApproximErrorNoisy(:,m), '.r-');
        legend(pl, 'Refrence noise free','Reference noisy');
        ylabel('Wave Approx. Error, dB');
        xlabel('Test number')
        grid on;

        if (AirArguments == 1)
            TitleStrings={'alfa0', 'freq0', 'c2', 'n', 'ro2', 'h', 'T', 'P'};
        else
            TitleStrings={'alfa0', 'freq0', 'c2', 'n', 'ro2', 'h'};
        end

        figure(figure_num); figure_num = figure_num + 1;
        for j = 1 : length(TitleStrings)
            subplot(length(TitleStrings) / 2, 2, j)
            plot(Args_Errors(:, j),'.b-');
            ylabel('Error, [%]');
            xlabel('Test num')
            title(['Error of ' TitleStrings{j}]);
            grid
        end % for j=1:length(TitleStrings)

    end % if (DrawFigures1==1)

end  %for m=1:Num_Parame_Error_Levels

ParNo = 2; % c1 parameter number 2
ArgNo = 6; % 6 - h parameter
%stdMedian_c1_Errors(m,:)=std(MedianArgs_Errors(:,ArgNo));






