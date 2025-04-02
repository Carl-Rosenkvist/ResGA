#!/bin/bash
#SBATCH --job-name=genfit
#SBATCH --partition=general1
#SBATCH --nodes=3
#SBATCH --ntasks=240
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=512
#SBATCH --time=07:00:00
#SBATCH --no-requeue
#SBATCH --mail-type=FAIL

# Set your desired parameters
threads=8
max_generations=30
population_size=16
shave=0.5

wd="/Users/carl/Phd/ResGA"
cross_dir="$wd/experimental_data/cross_sections"
emulator_path="$wd/Emulator-fast/build"
score_path="$wd/generations/scores.csv"
case="default"

# Run the Python script with the specified parameters
python3 ../pythonScripts/main.py --threads $threads --max_generations $max_generations --population_size $population_size --shave $shave --wd "$wd" --score_path "$score_path" --case $case

