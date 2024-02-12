function error = Objective_Func(args, ref_signal, meas_signal, params, f_axis)

model_signal = Wave_synthesize(args, ref_signal, params, f_axis);

% error = sum((model_signal(start_point : end_point) - meas_signal(start_point : end_point)) .^ 2) / ...
%         sum(meas_signal(start_point : end_point) .^ 2);

error = sum((model_signal - meas_signal) .^ 2) / sum(meas_signal .^ 2);

end


% Objective_Function = @(func_args)  Objective_Func(func_args, Wave_Ref, Wave_Meas, ...
%     acq_parameters_shifted, freq_axis, time_axis, ...
%     H_electronics, np1, np2, TRANSFER_MODEL); % model description function
% 
% [FittedArgs, approximation_error(i), exitflagPSO,outputPSO] = ...
%     particleswarm(Objective_Function, nvariables, LB, UB, options);