%% SETUP
clear; close all;

DrawFigures = 1; % 0 - Off; 1 - On; 2 - Draw each iteration results
figure_num = 1;

SaveToFile = 0;
savefile = ''; % ''
loadfile = ''; % ''

% 0 - use c1 and ro1
% 1 - use Temperature and Pressure
model_type = 0;

%% REFRENCE
% Generates a refrence chirp signal
N_points = 1024 * 2;
F_sampling = 10e6;
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

if (DrawFigures == 1)
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

if (model_type == 0)
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

if (model_type == 1)
    %% model_type == 1; use Temperature and Pressure model
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
end
Num_Parame_Error_Levels = length(e_c1);




