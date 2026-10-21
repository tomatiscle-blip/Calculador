@echo off
REM ============================================================
REM  Calculador - Ventana de la aplicacion (PySide6)
REM  Doble clic para abrir.
REM ============================================================
chcp 65001 >nul
cd /d "%~dp0"

if not exist ".venv\Scripts\python.exe" (
    echo No se encontro el entorno virtual en .venv
    pause
    exit /b 1
)

start "" pyw -m app

echo.
echo (La aplicacion se cerro. Codigo: %ERRORLEVEL%)
pause

