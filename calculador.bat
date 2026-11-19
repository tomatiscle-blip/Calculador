@echo off
REM ============================================================
REM  Calculador - Ventana de la aplicacion (PySide6)
REM  Doble clic para abrir.
REM ============================================================
chcp 65001 >nul
cd /d "%~dp0"

if not exist ".venv\Scripts\pythonw.exe" (
    echo No se encontro el entorno virtual en .venv
    pause
    exit /b 1
)

start "" "%~dp0.venv\Scripts\pythonw.exe" -m app
exit /b
