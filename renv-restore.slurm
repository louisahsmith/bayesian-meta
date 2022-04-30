#!/bin/bash
#SBATCH -n 1
#SBATCH -t 100 # Runtime in minutes
#SBATCH -p short # Partition to submit to
#SBATCH --mem=4000 # Memory per cpu in MB (see also --mem-per-cpu)
#SBATCH --open-mode=append
#SBATCH -o renv-restore.out # Standard out goes to this file
#SBATCH -e renv-restore.err # Standard err goes to this filehostname
#SBATCH --job-name=renv-restore

cd bayesian-meta
R CMD BATCH --quiet --no-restore --no-save R/renv-restore.R ../renv-restore.out
