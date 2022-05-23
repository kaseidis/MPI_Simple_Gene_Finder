#ifndef _GENE_JUDGE_H
#define _GENE_JUDGE_H

#include "./lib/Sequence.h"
#include "./lib/GeneRange.h"

#ifdef _MSC_VER // For MSVC
    #define CROSS_PLATFORM_HIDDEN_API
    #ifdef CROSS_PLATFORM_LIBRARY_EXPORTS
        #define CROSS_PLATFORM_API __declspec(dllexport)
    #else
        #define CROSS_PLATFORM_API __declspec(dllimport)
    #endif
#else // For GCC
    #define CROSS_PLATFORM_API __attribute((visibility("default")))
    #define CROSS_PLATFORM_HIDDEN_API __attribute((visibility("hidden")))
#endif

/**
 * @brief Check if a range in section of sequence is a gene.
 *        User can following this api do define their own gene finding
 *        criteria.
 * 
 * @param range 
 * @param seq  
 */
CROSS_PLATFORM_API bool isGene(const gene::GeneRange & range, const Sequence & seq);

#endif //_GENE_JUDGE_H