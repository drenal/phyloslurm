#!/bin/bash
# -------------------
# Slurm JOB ARRAY starter script fro phylobayes_starter.py
# -------------------
#SBATCH --job-name=start_phylobayes
#SBATCH --output=start_phylobayes_%A_%a.out
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
#
## JOB ARRAY settings
# example sbatch --array=1-100%50
# 1   ... starting job array subjob id
# 100 ... ending job array subjob id
# 50  ... number of concurrent job array subjobs
#SBATCH --array=1-100%50

# We assume that a directory with only alignment files was
# given to this script
# e.g. sbatch start_phylobayes.sbatch.sh /path/to/input_files_dir

EXTENSION=${2:-phylip}

FILE_TO_USE=`ls $1 | grep "${EXTENSION}$" | sed -n "${SLURM_ARRAY_TASK_ID}p"`

./start_phylobayes.py --input="${1}/${FILE_TO_USE}" \
--pb-args="-cat -gtr"            \
--chains=2                       \
--max-diff-bpcomp=0.3            \
--effective-size-tracecomp=50.0  \
--relative-diff-tracecomp=0.3    \
--max-sample-size=100000         \
--min-sample-size=6667           \
--sleep-range-min=3600           \
--sleep-range-max=3600           \
#--processes=1  # processes is infered from Slurm environemnt variables

