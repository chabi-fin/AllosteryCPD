#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --mem=500
#SBATCH --time=2:00:00
#SBATCH --partition=agkeller
#SBATCH --qos=standard
#SBATCH --job-name=cat_drive

module add GROMACS/2021.5-foss-2021b-CUDA-11.4.1-PLUMED-2.8.0
MDP="/home/lf1071fu/tcda/closed_holo/MDP"

files=$(ls -v run*/md_run/fitted_traj.xtc)
commands=$(printf 'c\n%.0s' {1..3})'\nc'

echo $files[@]
echo $commands
echo "Number of sims in concatenation ${#files[@]}"

# Concatenate trajectories
echo -e "$commands" | gmx trjcat -f ${files[@]} -o cat_trj.xtc -nobackup -settime

echo "24 24 24" | gmx trjconv -f cat_trj.xtc -s run1/md_run/closed_holo.tpr -center yes -pbc cluster -o centered_traj.xtc -n index.ndx -nobackup

echo "4 24" | gmx trjconv -f centered_traj.xtc -s run1/md_run/closed_holo.tpr -fit rot+trans -o fitted_traj.xtc -n index.ndx -nobackup

echo "24" | gmx trjconv -f fitted_traj.xtc -s run1/md_run/closed_holo.tpr -o fitted_traj_100.xtc -skip 100 -n index.ndx -nobackup

rm cat_trj.xtc centered_traj.xtc

