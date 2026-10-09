#!/bin/bash
# Build EasyCGTree_GUI as a program for Linux (run on Linux, Python 3.9 or newer).
# The result is 'dist/EasyCGTree_GUI'; copy it into the EasyCGTree directory (next to EasyCGTree.pl).
set -e
cd "$(dirname "$0")"
python3 -m venv .venv
. .venv/bin/activate
pip install --upgrade pip
pip install "PySide6-Essentials>=6.5" "pyinstaller>=6"
pyinstaller --noconfirm --clean --onefile --windowed --name EasyCGTree_GUI ../EasyCGTree_GUI.py
echo
echo "Done: $(pwd)/dist/EasyCGTree_GUI"
