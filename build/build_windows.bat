@echo off
rem Build EasyCGTree_GUI.exe for Windows (run on Windows, Python 3.9 or newer from python.org).
rem The result is 'dist\EasyCGTree_GUI.exe'; copy it into the EasyCGTree directory (next to EasyCGTree.pl).
cd /d "%~dp0"
py -3 -m venv .venv || goto :error
call .venv\Scripts\activate.bat || goto :error
python -m pip install --upgrade pip || goto :error
pip install "PySide6-Essentials>=6.5" "pyinstaller>=6" || goto :error
pyinstaller --noconfirm --clean --onefile --windowed --name EasyCGTree_GUI ..\EasyCGTree_GUI.py || goto :error
echo.
echo Done: %cd%\dist\EasyCGTree_GUI.exe
pause
exit /b 0
:error
echo Build failed.
pause
exit /b 1
