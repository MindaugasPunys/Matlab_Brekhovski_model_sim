function ret = plot_xilinx_data(dir, var)
% Construct the file path
file_path = dir;

% Read the file
fid = fopen(file_path, 'r');
if fid == -1
    error('Could not open the file');
end

% Initialize the data variable
data = [];

% Read the file line by line
line = fgetl(fid);
while ischar(line)
    % Check if the line contains the variable
    if contains(line, var)
        % Extract the data from the line
        % Assuming the format is "var[index]: value1, value2"
        var_data = sscanf(line, [var '[%d]: %f, %f']);
        % Append the data to the result
        data = [data; var_data(2:end)']; % Exclude the index
    end
    line = fgetl(fid);
end

% Close the file
fclose(fid);

% Check if data is empty
if isempty(data)
    error('Variable %s not found in the log file.', var);
end

% Plot the data
figure;
plot(data);
title(sprintf('Plot of %s', var));
xlabel('Index'); ylabel('Value');
grid on;

ret = data;
end
