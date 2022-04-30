#!/bin/bash
#SBATCH -n 8
#SBATCH -t 20 # Runtime in minutes
#SBATCH -p debug # Partition to submit to
#SBATCH --mem=4000 # Memory per cpu in MB (see also --mem-per-cpu)
#SBATCH --open-mode=append
#SBATCH -o test-stan.out # Standard out goes to this file
#SBATCH -e test-stan.err # Standard err goes to this filehostname
#SBATCH --job-name=test-stan

cd bayesian-meta
R CMD BATCH --quiet --no-restore --no-save R/run-all.R ../test-stan.out
