cmake_minimum_required(VERSION 3.13)

project(GeneFinder)
#include(CheckLanguage)
#check_language(CUDA)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

find_package(OpenMP)
find_package(MPI)

# Gene gene_judge library
add_library(gene_judge SHARED ./gene_judge/gene_judge.cpp ./gene_judge/lib/Sequence.cpp)
target_compile_features(gene_judge PRIVATE cxx_std_17)

# Non MPI Version
add_executable(gene_finder ./src/main.cpp ./src/lib/orf_finder.cpp ./src/lib/Sequence.cpp ./src/lib/Fasta.cpp ./src/lib/InputParser.cpp)
target_link_libraries (gene_finder gene_judge)
if (OPENMP_FOUND)
    if (NOT WIN32)
        target_link_libraries(gene_finder OpenMP::OpenMP_CXX m)
    else()
        target_link_libraries(gene_finder OpenMP::OpenMP_CXX)
    endif()
endif()
target_compile_features(gene_finder PRIVATE cxx_std_17)

# MPI Version
if (MPI_FOUND)
    add_executable(gene_finder_mpi ./src/main_mpi.cpp ./src/lib/orf_finder.cpp ./src/lib/Sequence.cpp ./src/lib/Fasta.cpp ./src/lib/InputParser.cpp)
    include_directories(SYSTEM ${MPI_INCLUDE_PATH})
    target_link_libraries (gene_finder_mpi gene_judge)
    target_link_libraries(gene_finder_mpi ${MPI_CXX_LIBRARIES})
    if (OPENMP_FOUND)
        if (NOT WIN32)
            target_link_libraries(gene_finder_mpi OpenMP::OpenMP_CXX m)
        else()
            target_link_libraries(gene_finder_mpi OpenMP::OpenMP_CXX)
        endif()
    endif()
    target_compile_features(gene_finder_mpi PRIVATE cxx_std_17)
endif()

#if (CMAKE_CUDA_COMPILER)
#    enable_language(CUDA)
#    add_executable(ray_trace_cuda ray_trace.cu bitmap.c timer.c)
#    if (NOT WIN32)
#        target_link_libraries(ray_trace_cuda m)
#    endif()
#endif ()