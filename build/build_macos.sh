#!/bin/bash
# Build EasyCGTree_GUI as a macOS application (run on a Mac, Python 3.9 or newer, e.g. from python.org).
# The result is 'dist/EasyCGTree_GUI.app'; copy it into the EasyCGTree directory (next to EasyCGTree.pl).
# A Mac with Apple silicon (M1...) builds an app for Apple silicon; an Intel Mac builds one for Intel.
set -e
cd "$(dirname "$0")"
python3 -m venv .venv
. .venv/bin/activate
pip install --upgrade pip
pip install "PySide6-Essentials>=6.5" "pyinstaller>=6"
pyinstaller --noconfirm --clean --windowed --name EasyCGTree_GUI ../EasyCGTree_GUI.py
echo
echo "Done: $(pwd)/dist/EasyCGTree_GUI.app"
echo "The app is not signed: at the first start, right-click it and choose 'Open'."
