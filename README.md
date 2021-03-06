# MPI Simple Gene Finder Framework

# Summary

In bioinformatics, gene finding refers to identifying regions of genomic DNA that encode genes. Gene finding includes protein-coding genes and RNA genes. Gene discovery is the first and one of the most critical steps in understanding the sequencing of a species' genome. However, the process of gene finding takes a lot of time. This program is a framework that can help us do gene finding job on multiple nodes using MPI.

## Compile

### Pre-request:

- GCC/Clang++ (Support C++17/OpenMP)
- OpenMPI Development pack
- CMake
- Make

For Ubuntu:
```
sudo apt install cmake make g++ openmpi-bin libopenmpi-dev
```

### Compile Command
```
cmake .
make
```

## Built Binary
There is a X64 Linux Build (openMPI) in folder: [./build/](./build/)

## Run
Single Node Version:
```
Usage: ./gene_finder --input INPUT_FILE_PATH --output OUTPUT_FILE_PATH [--pattern LABEL_PATTERN --output-line-width WIDTH --time]
    Default:
        LABEL_PATTERN = '%s | gene | frame=%d | LOC=[%d,%d]'
        WIDTH = 70
```

Mutiple Node (MPI) Versoin:
```
Usage: mpirun [MPI_ARGS] ./gene_finder_mpi --input INPUT_FILE_PATH --output OUTPUT_FILE_PATH [--pattern LABEL_PATTERN --output-line-width WIDTH]
    Default:
        LABEL_PATTERN = '%s | gene | LOC=[%d,%d]'
        WIDTH = 70
```

Here are sample run command sbatch script:
- [Single Node Version](./build/run_gene_finder.sh)
- [MPI Version](./build/run_gene_finder_mpi.sbatch)


## Development

### Self Defined ORF Evaluation Function

All code for developing library are in [``./gene_judge/``](./gene_judge/) directory.

You can edit [``./gene_judge/gene_judge.cpp``](./gene_judge/gene_judge.cpp) to program you own gene identifing function.

Here are two other sample:

- [``gene_judge_filter_all.cpp``](./gene_judge/gene_judge_filter_all.cpp): Filter out all orfs.
- [``gene_judge_get_all.cpp``](./gene_judge/gene_judge_get_all.cpp): Get all orfs.

### Compile Self Defined ORF Evaluation Function
Go root directory of this repo, then:
```
cmake .
make gene_judge
```
It will genreate a dynamic linked library file. You can replace the file in build directory (.so or .dll) to the one you build.

## Paper & Presntation

[``Distributed Framework for Gene Finding using Open-MPI``](./paper/paper.pdf)

[``Presentation Slides``](./paper/presentation.pdf)