#!/bin/bash
cat modules/InputParser.py > chemkin.py
echo '' >> chemkin.py
echo '' >> chemkin.py
cat modules/ReactionCoeffs.py >> chemkin.py
echo '' >> chemkin.py
echo '' >> chemkin.py
cat modules/chemkin.py >> chemkin.py