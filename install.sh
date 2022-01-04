#!/bin/bash

#python3 -m PyInstaller --onefile eval.py
pyinstaller --onefile --hidden-import [numpy,numpy.distutils] eval.py
#pyinstaller --additional-hooks-dir=hooks --onefile eval.py
