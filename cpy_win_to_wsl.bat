@echo off
set "source=D:\Magistras\Magistro projektas\Matlab_realization\Raw\test"
set "destination=\\wsl.localhost\Ubuntu\home\user\projects\pso_model\test"

echo Copying files from %source% to %destination% ...

rem Check if source directory exists
if not exist "%source%" (
    echo Source directory does not exist.
    pause
    exit /b 1
)

rem Create destination directory if it doesn't exist
mkdir "%destination%" 2>nul

rem Use robocopy to copy files
robocopy "%source%" "%destination%" /e

rem Display success message
echo Files copied successfully.
REM pause
