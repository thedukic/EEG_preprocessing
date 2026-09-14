@echo off
setlocal enabledelayedexpansion

:: Define the target directory
set "target_dir=E:\1_EEG_DATA\ALS\T1"

:: Check if the target directory exists
if not exist "%target_dir%" (
    echo Target directory %target_dir% does not exist.
    exit /b
)

echo Checking subfolders in %target_dir%...
echo ===================================================

:: Loop through each subfolder
for /d %%D in ("%target_dir%\*") do (
    set "missing="
    
    :: Check for each required suffix
    for %%S in (_EC1 _EC2 _EC3 _EO1 _EO2 _EO3) do (
        dir /b /a-d "%%D\*%%S*" >nul 2>&1
        if errorlevel 1 (
            set "missing=1"
        )
    )
    
    :: If any of the files were missing, print the subfolder name
    if defined missing (
        echo Missing required files in: %%~nxD
    )
)

echo ===================================================
echo Check complete.
pause