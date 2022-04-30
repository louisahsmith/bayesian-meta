#!/bin/bash
#SBATCH -n 64
#SBATCH -t 400 # Runtime in minutes
#SBATCH -p short # Partition to submit to
#SBATCH --mem=4000 # Memory per cpu in MB (see also --mem-per-cpu)
#SBATCH --open-mode=append
#SBATCH -o Historical04292022.out # Standard out goes to this file
#SBATCH -e Historical04292022.err # Standard err goes to this filehostname
#SBATCH --job-name=Historical04292022

cd bayesian-meta
R CMD BATCH --quiet --no-restore --no-save R/run-all.R ../bayesian-meta.out
