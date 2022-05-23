#pragma once
#ifndef _ORF_FINDER_H
#define _ORF_FINDER_H
#include <vector>
#include <string>
#include <stdio.h>
#include "Fasta.h"
#ifdef __GNUC__
#define PACK( __Declaration__ ) __Declaration__ __attribute__((__packed__))
#endif

#ifdef _MSC_VER
#define PACK( __Declaration__ ) __pragma( pack(push, 1) ) __Declaration__ __pragma( pack(pop))
#endif
namespace gene
{

     /**
      * @brief Range of certain sequence
      */
     typedef struct _gene_range
     {
          // using unsigned long long for openMPI
          /**
           * @brief start postion of gene
           */
          unsigned long long start;
          /**
           * @brief end position of gene
           */
          unsigned long long end;
          /**
           * @brief frame position of gene
           */
          int8_t frame;
          /**
           * @brief Get length of gene
           * 
           * @return size_t 
           */
          inline unsigned long long length() const
          {
               return this->abs_end() - this->abs_start() + 1;
          }
          /**
           * @brief Get real start position of gene range (smaller one)
           * 
           * @return size_t 
           */
          inline unsigned long long abs_start() const
          {
               return frame < 0 ? this->end : this->start;
          }
          /**
           * @brief Get real end position of gene range (larger one)
           * 
           * @return size_t 
           */
          inline unsigned long long abs_end() const
          {
               return frame < 0 ? this->start : this->end;
          }
     } GeneRange;
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