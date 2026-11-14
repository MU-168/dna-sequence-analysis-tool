@echo off
rem DNA analysis teaching tool - GUI launcher
rem Double-click this file to open the input window.
setlocal
cd /d "%~dp0"

where pyw >nul 2>nul
if errorlevel 1 goto try_py
start "" pyw -3 "dna_sequence_analysis.py"
exit /b 0

:try_py
where py >nul 2>nul
if errorlevel 1 goto no_python
start "" py -3 "dna_sequence_analysis.py"
exit /b 0

:no_python
echo Python launcher not found. Install Python 3 and run: pip install biopython
pause
