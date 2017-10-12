#!/bin/bash
cat InputParser.py > chemkin.py
echo '' >> chemkin.py
cat reaction_coeffs.py >> chemkin.py
echo '' >> chemkin.py
cat reaction_rates.py >> chemkin.py