#include "gene_judge.h"
#include <omp.h>

/**
 * Judge ORF is a gene or not, this is a simple demo that following this
 * criteria:
 *  - Must begin from a start codon
 *  - ORF must contain at least 96 bp (32 amino acids)
 *  - ORF must occur in a CpG
 *  - CpG should appear in first third of the sequence
 * 
 * User can modifie this file, then using ``make gene_judge`` command
 * that process gene finding algorithm by they defined.
 */
bool isGene(const gene::GeneRange &range, const Sequence &seq)
{
    auto start = range.abs_start();
    auto end = range.abs_end();
    end = start + range.length() / 3;
    start = end-200 < start? end-200: start;
    start = start > 0? start:0;

    int n = 200; double t_ratio = 0.6; double t_gc = 0.5;

    std::string_view seq_view(
        seq.getSequence().c_str(),
        seq.getSequence().length());
    auto l = seq_view.length();
    if (l < n || start >= l - n || end > l - n)
    {
        return false;
    }
    bool result = false;
    // Searching Cpg island
    #pragma omp parallel for
    for (int64_t i = start; i < end; ++i)
    {
        // Getting nC, nG, nCpG for current window
        size_t nC = 0, nG = 0, nCpG = 0;
        for (int64_t j = i; j < i + 200; ++j)
        {
            if (seq_view[j] == 'C')
            {
                ++nC;
                if (j + 1 < l && seq_view[j + 1] == 'G')
                    ++nCpG;
            }
            else if (seq_view[j] == 'G')
                ++nG;
        }
        // Get Obs/Exp and GC content
        double oe_ratio = nCpG;
        oe_ratio = oe_ratio / (nC * nG) * n;
        double gc_content = nC + nG;
        gc_content /= n;
        // Check and set result, not return because OpenMP
        if (oe_ratio > t_ratio && gc_content > t_gc)
            result = true;
    }

    return result;
}