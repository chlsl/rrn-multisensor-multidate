#!/usr/bin/env bash

normalization=$(dirname "$0")/../rrn/normalization.py
data=$(ls $(dirname "$0")/data/versailles-short/{S2/L1C,S2/L2A,L8/pansharp-upsa}/*band_B0{2,3,4}*.tif)

echo ""
echo "=== Applying $normalization with three bands: B02, B03 and B04"
echo ""

python3 $normalization $data --band B04 B03 B02 --outdir out

