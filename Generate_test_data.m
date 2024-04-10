%% SETUP
clear; close all;

DrawFigures = 1; % 0 - Off; 1 - On; 2 - Draw each iteration results
figure_num = 1;

SaveToFile = 0;
savefile = ''; % ''
loadfile = ''; % ''

FixedNoise = 0; % load wave record from loadfile

% 0 - use c1 and ro1
% 1 - use Temperature and Pressure (Better results)
global AirArguments;
AirArguments = 1;
TestCycles_New = 5;  % times to repeat aproximation

%% REFRENCE
% Generates a refrence chirp signal
N_points = 1024; % 2048
F_sampling = 4e6; % 10e6
T_sampling = 1 / F_sampling;
Magnitude = 1000;

% https://en.wikipedia.org/wiki/Speed_of_sound
% Temperature = 20; % [degC]
% Pressure = 101.325e3; % [Pa] standard atmospheris pressure
% c1 = 331.38 * sqrt(1 + Temperature / 273.15); % speed of sound
% gama = 1.4; % Air adiabatic constant
% Ks = gama * Pressure; % Coefficient of stiffness
% ro1 = Ks / (c1^2);

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

if (DrawFigures >= 1)
    %% Plot the generated signal
    figure(figure_num); figure_num = figure_num + 1;

    plot(time, Reference);
    xlabel('Time (s)');
    ylabel('Amplitude');
    title('Chirp Excitation Signal');
end


%% MODEL PARAMETERS

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

    OriginalArgs = [alfa0, freq0, V_sluoksnio, n, ro2, h];
    LB_coeff = -[50, 99, 50, 50, 50, 99]; % [%]
    UB_coeff = [50, 100, 50, 50, 50, 100]; % [%]
    acquisition_parameters = [F_sampling, c1, ro1];
    acquisition_parameters_error = 0; % [%]
    e_c1 = 0;
else
    %% AirArguments == 0; use c1 and ro1 model

    c1 = 343; % First media ultrasound velocity (Air)
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

if (DrawFigures >= 1)
    %% Plot the generated ultrasound wave
    figure(figure_num); figure_num = figure_num + 1;
    plot(Wave_Meas_WithoutNoise,'r');
    title('Refrence signal without noise');
    grid
end

% GaindB = 60;  % [dB]
% Gain = 10^(GaindB/20);
Gain = max(Reference) / max(Wave_Meas_WithoutNoise);


%% EXPORT DATA

T = TransferFunction(OriginalArgs, acquisition_parameters, freq_axis);
Wave = Wave_synthesize(OriginalArgs, Reference, acquisition_parameters, freq_axis);
% IN
export_data(OriginalArgs, 'test/TransferFunction/OriginalArgs');
export_data(acquisition_parameters, 'test/TransferFunction/acquisition_parameters');
export_data(freq_axis, 'test/TransferFunction/freq_axis');
export_data(Reference, 'test/TransferFunction/Reference');
export_data(Wave_Meas_WithoutNoise, 'test/TransferFunction/Measured');
% OUT
export_data(T, 'test/TransferFunction/T');
export_data(Wave, 'test/TransferFunction/Wave');

% COPY TEST DATA TO WSL
system('cpy_win_to_wsl.bat');

%% TEST

% dir = 'C:/GIT/VITIS/Model_sim/solution/csim/report/TransferFunction_csim.log';
% dir = 'C:/GIT/VITIS/Model_sim/solution/csim/report/WaveSynthesis_csim.log';
% dir = 'C:/GIT/VITIS/Model_sim/solution/csim/report/pso_process_csim.log';
dir = '//wsl.localhost/Ubuntu/home/user/projects/pso_model/log.txt';

% TransferFunction
plot_xilinx_data(dir, 'T');
hold on;
plot(real(T), '--')
hold on;
plot(imag(T), '--')
title('T koeficientai'); legend('HLS', 'Matlab');
xlabel('Indeksas'); ylabel('Amplitudė');

%FFT wrapper module
FFT_TEST = 0
if FFT_TEST
    clear; close all;
    dir_fft = 'C:/GIT/VITIS/FFT_TEST/solution1/csim/report/FFT_call_csim.log';

    fft_input = plot_xilinx_data(dir_fft, 'fft_input'); xlim([0 100]);
    
    fft_output = plot_xilinx_data(dir_fft, 'fft_output');
    matlab_fft = fft(fft_input);
    hold on; plot(real(matlab_fft), '--');
    hold on; plot(imag(matlab_fft), '--');
    legend('HLS re','HLS im','M re','M im')

    ifft_output = plot_xilinx_data(dir_fft, 'ifft_output'); xlim([0 100]);

    figure;
    plot(fft_input);
    hold on;
    plot(ifft_output, '--');
    legend('fft input', 'ifft output'); xlim([0 100]);
end

% WaveSynthesize
plot_xilinx_data(dir, 'wave_out');
hold on;
plot(Wave, '--')
title('Wave'); legend('HLS', 'Matlab');
xlabel('Indeksas'); ylabel('Amplitudė');

% PSO
hls_args = plot_xilinx_data(dir, 'hls_args');
hls_wave = Wave_synthesize(hls_args, Reference, acquisition_parameters, freq_axis);
figure;
plot(hls_wave, '-')
hold on;
plot(Wave_Meas_WithoutNoise, '--')

%%
% HDL coder result?
