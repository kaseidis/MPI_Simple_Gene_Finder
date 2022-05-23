#!/bin/bash
#SBATCH --nodes=1
#SBATCH --time=0:01:00
#SBATCH --job-name=gene_finder
#SBATCH --partition=fast
#SBATCH --error=./gene_finder.err
#SBATCH --output=./gene_finder.out
#SBATCH --ntasks-per-node=1
#SBATCH --exclusive

./gene_finder --input ./data/chr21.fasta --output ./chr21.fasta.out --time