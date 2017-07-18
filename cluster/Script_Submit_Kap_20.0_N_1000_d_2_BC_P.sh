#!/bin/bash
#PBS -N RGG_1000_20.0_2_P
#PBS -m be
#PBS -q standard
cd /scratchcomp03/MGarrod_4/
chmod u+x /scratchcomp03/MGarrod_4/RGG_Sample_Mu2_Only.py
python /scratchcomp03/MGarrod_4/RGG_Sample_Mu2_Only.py 1000 20.0 2 P