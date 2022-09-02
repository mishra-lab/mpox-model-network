#!/bin/bash
#SBATCH --job-name mpxv-netcb
#SBATCH --nodes=1
#SBATCH --cpus-per-task=40
#SBATCH --time=1:00:00

module load gcc/9.2.0
module load r/3.6.3
# install.packages( # might be missing some
#   c('ggplot2','ggExtra','ggridges','viridis','reshape2'))
Rscript scinet/main.r

# using this script on scinet:
# sbatch submit.sh
