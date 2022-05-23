#pragma once
#ifndef _ORF_FINDER_H
#define _ORF_FINDER_H
#include <vector>
#include <string>
#include <stdio.h>
#include "Fasta.h"
#include "GeneRange.h"
namespace gene
{    
     /**
      * @brief Get orfs from dna/rna sequence, returns vector of GeneRange object.
      * 
      * @param seq 
      * @param frame 
      * @param startLoc 
      * @param endLoc 
      * @return std::vector<GeneRange> 
      */
     std::vector<GeneRange> getORFS(
         const Sequence &seq, int8_t frame, size_t startLoc,
         size_t endLoc);
}
#endif