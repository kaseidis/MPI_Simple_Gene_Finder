#!/bin/bash
# A sample sbatch file for gene_finder_mpi


#SBATCH --job-name=gene_finder
#SBATCH --output=gene_finder_mpi.out
#SBATCH --ntasks=5
#SBATCH --partition=general

# Load the default OpenMPI module.
# module load openmpi

mpirun ./gene_finder_mpi --input ./data/chr21.fasta --output ./chr21_mpi.fasta.out