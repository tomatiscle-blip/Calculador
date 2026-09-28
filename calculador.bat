@echo off
REM ============================================================
REM  Calculador - Ventana de la aplicacion (PySide6)
REM  Doble clic para abrir.
REM ============================================================
chcp 65001 >nul
cd /d "%~dp0"
start "" pyw -m app
if errorlevel 1 start "" py -m app
