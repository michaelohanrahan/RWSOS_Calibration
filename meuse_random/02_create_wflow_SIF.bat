@echo off
REM Activate WSL with Ubuntu distribution and change permissions of wflow_container.sh
wsl -d Ubuntu -e bash -c "chmod +x wflow_container.sh"

REM Execute the wflow_container.sh script
wsl -d Ubuntu -e bash -c "./wflow_container.sh"