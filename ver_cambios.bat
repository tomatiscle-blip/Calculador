@echo off
REM ============================================================
REM  Que cambio en el proyecto - doble clic
REM  Lista los archivos nuevos y los modificados, en castellano.
REM ============================================================
chcp 65001 >nul
cd /d "%~dp0"
py tools\ver_cambios.py
echo.
echo ------------------------------------------------------------
echo  Para abrir los archivos nuevos en VS Code:
echo      py tools\ver_cambios.py --abrir
echo  Para comparar un archivo antes/despues:
echo      py tools\ver_cambios.py --diff P06_Portico_dxf.py
echo ------------------------------------------------------------
pause
