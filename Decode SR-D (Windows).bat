@echo off
cd /d "%~dp0"
python -c "import numpy, reedsolo, cv2, rawpy" 2>nul || python -m pip install --user -r requirements.txt
start "" pythonw srd_gui.py
