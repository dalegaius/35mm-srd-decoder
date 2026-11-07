#!/bin/bash
cd "$(dirname "$0")"
python3 -c "import numpy, reedsolo, cv2, rawpy" 2>/dev/null || python3 -m pip install --user -r requirements.txt
python3 srd_gui.py
