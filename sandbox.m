clc; clear;
% Generate random input data
n = 1000; % Number of data points
model_signal = randn(n, 1); % Random model signal
meas_signal = randn(n, 1); % Random measurement signal

% Approach 1: Original MATLAB calculation
error_matlab = sum((model_signal - meas_signal) .^ 2) / sum(meas_signal .^ 2);

% Approach 2: Custom C-like calculation
numerator = 0;
denominator = 0;
for i = 1:n
    diff = model_signal(i) - meas_signal(i);
    numerator = numerator + diff^2;
    denominator = denominator + meas_signal(i)^2;
end
error_c_like = numerator / denominator;

% Compare the results
tolerance = 1e-10; % Tolerance for floating-point comparison
if abs(error_matlab - error_c_like) < tolerance
    disp('Results are comparable.');
else
    disp('Results are not comparable.');
end
