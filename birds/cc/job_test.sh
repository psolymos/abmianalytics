#!/bin/bash
#SBATCH --account=def-psolymos
#SBATCH --ntasks=8
#SBATCH --mem-per-cpu=2048M
#SBATCH --time=00:30:00
#SBATCH --job-name=test_makecluster
#SBATCH --output=%x-%j.out
#SBATCH --mail-user=solymos@ualberta.ca
#SBATCH --mail-type=ALL

# Load modules
module nixpkgs/16.09
module load gcc/7.3.0
module load openmpi/3.1.2
module load r/3.5.1

# Export the nodes names.
# If all processes are allocated on the same node,
# NODESLIST contains : node1 node1 node1 node1
export NODESLIST=$(echo $(srun hostname))

# Run R script
Rscript --vanilla test.R
