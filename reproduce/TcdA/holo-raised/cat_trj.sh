#!/bin/bash
#SBATCH --time=7-00:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=20G
#SBATCH --partition=normal
#SBATCH --job-name=cat_drive

module add bio/GROMACS/2021.5-foss-2021b-CUDA-11.4.1-PLUMED-2.8.0

MDP="/scratch/hpc-prf-cpdallo/tcda/open_holo/MDP"

files=$(ls -v run*/md_run/fitted_traj.xtc)
commands=$(printf 'c\n%.0s' {1..3})'\nc'

echo $files[@]
echo $commands
echo "Number of sims in concatenation ${#files[@]}"

# Concatenate trajectories
echo -e "$commands" | gmx trjcat -f ${files[@]} -o cat_trj.xtc -nobackup -settime

echo "24 24 24" | gmx trjconv -f cat_trj.xtc -s run1/md_run/open_holo.tpr -center yes -pbc cluster -o centered_traj.xtc -n index.ndx -nobackup

echo "4 24" | gmx trjconv -f centered_traj.xtc -s run1/md_run/open_holo.tpr -fit rot+trans -o fitted_traj.xtc -n index.ndx -nobackup

echo "24" | gmx trjconv -f fitted_traj.xtc -s run1/md_run/open_holo.tpr -o fitted_traj_100.xtc -skip 100 -n index.ndx -nobackup

rm cat_trj.xtc centered_traj.xtc

