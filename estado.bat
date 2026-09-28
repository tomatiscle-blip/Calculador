@echo off
REM ============================================================
REM  Estado del proyecto - doble clic para ver el semaforo
REM  (calculos al dia / faltantes / a desarrollar)
REM ============================================================
chcp 65001 >nul
cd /d "%~dp0"
py -m calc.pipeline %*
echo.
pause
