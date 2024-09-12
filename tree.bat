@echo off
setlocal enabledelayedexpansion

set "outputFile=folder_tree.txt"
set "targetLevel=2"
set "currentLevel=0"

(
    call :printTree "." 0
) > "%outputFile%"
exit /b

:printTree
set "folderPath=%~1"
set "level=%~2"

for /d %%D in ("%folderPath%\*") do (
    set "folderName=%%~nxD"
    set "indent="
    for /L %%I in (1,1,!level!) do set "indent=!indent!    "
    echo !indent!!folderName!
    
    if !level! lss !targetLevel! (
        call :printTree "%%D" !level!+1
    )
)
exit /b