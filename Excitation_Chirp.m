function [Reference, time] = Excitation_Chirp(F1_chirp, F2_chirp, T_chirp, T_StartShift, Magnitude, XducF1, XducF2, T_sampling, T_All)
    % Calculate time vector
    time = 0:T_sampling:T_All;

    % Shift time vector
    time = time + T_StartShift;

    % Generate chirp signal
    Chirp_signal = chirp(time, F1_chirp, T_chirp, F2_chirp);

    % Apply excitation factors
    Reference = Magnitude * (XducF1 + XducF2) .* Chirp_signal;
end
