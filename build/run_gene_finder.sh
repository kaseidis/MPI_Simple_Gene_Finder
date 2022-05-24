#!/bin/bash
#SBATCH --nodes=1
#SBATCH --time=0:01:00
#SBATCH --job-name=gene_finder
#SBATCH --partition=fast
#SBATCH --error=./output/slurm/gene_finder_%j.err
#SBATCH --output=./output/times/gene_finder_%j.out
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=2
#SBATCH --exclusive

export OMP_NUM_THREADS=2
export LD_LIBRARY_PATH=$(pwd):$LD_LIBRARY_PATH
./gene_finder --input ../data/chr21.fasta --output ./output/chr21.fasta.out --time