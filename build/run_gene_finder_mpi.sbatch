#!/bin/bash
#SBATCH --time=0:01:00
#SBATCH --job-name=gene_finder_mpi
#SBATCH --partition=fast
#SBATCH --error=./output/slurm/gene_finder_mpi_%j.err
#SBATCH --output=./output/times/gene_finder_mpi_%j.out
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=2
#SBATCH --exclusive

export OMP_NUM_THREADS=2
export LD_LIBRARY_PATH=$(pwd):$LD_LIBRARY_PATH
mpirun ./gene_finder_mpi --input ../data/chr21.fasta --output ./output/chr21_mpi.fasta.out --time