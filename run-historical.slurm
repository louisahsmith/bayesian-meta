#!/bin/bash
#SBATCH -n 64
#SBATCH -t 840 # Runtime in minutes
#SBATCH -p short # Partition to submit to
#SBATCH --mem=4000 # Memory per cpu in MB (see also --mem-per-cpu)
#SBATCH --open-mode=append
#SBATCH -o Historical012022.out # Standard out goes to this file
#SBATCH -e Historical012022.err # Standard err goes to this filehostname
#SBATCH --job-name=Historical012022

cd bayesian-meta
R CMD BATCH --quiet --no-restore --no-save R/eumaeus-historical.R ../bayesian-meta.out
