#pragma once
#ifndef _ORF_FINDER_H
#define _ORF_FINDER_H
#include <vector>
#include <string>
#include <stdio.h>
#include "Fasta.h"

namespace gene
{
    typedef struct _gene_range
    {
        // using unsigned long long for openMPI
        unsigned long long start;
        unsigned long long end;
        int8_t frame;
        inline size_t length() const{
             return this->abs_end() - this->abs_start() + 1;
        }
        inline size_t abs_start() const{
             return frame < 0?this->end:this->start;
        }
        inline size_t abs_end() const{
             return frame < 0?this->start:this->end;
        }
    } GeneRange;
    std::vector<GeneRange> getORFS(
        const Sequence &seq, int8_t frame, size_t startLoc, 
        size_t endLoc);
}
#endif