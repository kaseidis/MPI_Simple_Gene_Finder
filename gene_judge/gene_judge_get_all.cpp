#include "./lib/gene_judge.h"
#include <omp.h>
#include <string_view>

/**
 * Fillter all of orfs, always return invalid geneRange;
 */
gene::GeneRange isGene(const gene::GeneRange &range, const Sequence &seq)
{
    gene::GeneRange result{range.start,range.end,range.frame};
    return result;
}