#!/bin/bash
# -------------------
# Slurm starter script for start_phylobayes.py
# -------------------
#SBATCH --job-name=start_phylobayes
#SBATCH --output=start_phylobayes_%j.out
##SBATCH --partition=low
##SBATCH --mem-per-cpu=2048
##SBATCH --time=7-12:00:00
# Set ntasks to:
# number_of_processes [times] number_of_chains
# if number_of_processes == 1 -> pb
# if number_of_processes > 1 -> pb_mpi
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=1
##SBATCH --mail-type=END,FAIL
##SBATCH --mail-user=email@address.address

# We assume that a single (alignment) filename was
# given to this script
# e.g. sbatch start_phylobayes.sbatch.sh /path/to/input_file

./start_phylobayes.py --input=$1 \
--pb-args="-cat -gtr"            \
--chains=2                       \
--max-diff-bpcomp=0.3            \
--effective-size-tracecomp=50.0  \
--relative-diff-tracecomp=0.3    \
--max-sample-size=100000         \
--min-sample-size=6667           \
--sleep-range-min=3600           \
--sleep-range-max=3600           \
#--processes=1  # Processes is infered from Slurm environment variables

