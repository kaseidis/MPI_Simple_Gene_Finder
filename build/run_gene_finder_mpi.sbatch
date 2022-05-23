#!/bin/bash
#SBATCH --nodes=4
#SBATCH --time=0:01:00
#SBATCH --job-name=gene_finder_mpi
#SBATCH --partition=fast
#SBATCH --error=./gene_finder_mpi.err
#SBATCH --output=./gene_finder_mpi.out
#SBATCH --ntasks-per-node=1
#SBATCH --exclusive

mpirun ./gene_finder_mpi --input ../data/chr21.fasta --output ./chr21_mpi.fasta.out --time