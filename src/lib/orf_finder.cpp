#include "orf_finder.h"
#include <set>
#include <string_view>
#include <memory>
#include <string>
#include <stdexcept>

/**
 * @brief Check if a sequence is DNA sequence
 *
 * @param seq
 * @return true
 * @return false
 */
bool isDNA(const std::string &seq)
{
    auto l = seq.length();
    for (int i = 0; i < l; ++i)
        if (seq[i] == 'T')
            return true;
    return false;
}

/**
 * @brief Check if a sequence is RNA sequence
 *
 * @param seq
 * @return true
 * @return false
 */
bool toRNA(std::string &seq)
{
    auto l = seq.length();
    #pragma omp parallel for
    for (int i = 0; i < l; ++i)
        if (seq[i] == 'T')
            seq[i] = 'U';
    return true;
}

/**
 * @brief Convert RNA sequence to DNA sequence
 *
 * @param seq
 * @return true
 * @return false
 */
bool toDNA(std::string &seq)
{
    auto l = seq.length();
    #pragma omp parallel for
    for (int i = 0; i < l; ++i)
        if (seq[i] == 'U')
            seq[i] = 'T';
    return true;
}

/**
 * @brief Convert RNA sequence to 3'5 frame
 *
 * @param seq
 */
bool to35RNA(std::string &seq)
{
    auto l = seq.length();
    #pragma omp parallel for
    for (int i = 0; i < l; ++i)
        switch (seq[i])
        {
        case 'A':
            seq[i] = 'U';
            break;
        case 'T':
        case 'U':
            seq[i] = 'A';
            break;
        case 'C':
            seq[i] = 'G';
            break;
        case 'G':
            seq[i] = 'C';
            break;
        }
    return true;
}

/**
 * @brief Reverse a sequence
 *
 * @param seq
 */
void reverse(std::string &seq)
{
    auto l = seq.length();
    #pragma omp parallel for
    for (int i = 0; i < l / 2; i++)
    {
        char tmp = seq[i];
        seq[i] = seq[l - i - 1];
        seq[l - i - 1] = tmp;
    }
}

std::vector<gene::GeneRange> gene::getORFS(
    const Sequence &seq, int8_t frame, size_t startLoc,
    size_t endLoc)
{
    endLoc -= 1;
    // Get length
    const auto l = seq.getSequence().length();

    // Check for valid frame value
    if (frame == 0 || frame > 3 || frame < -3)
    {
        throw 1;
    }
    // Get duplicate sequence data
    std::string seqData = seq.getSequence();
    // Convert DNA to RNA
    bool dnaFlag = isDNA(seqData);
    if (dnaFlag)
        toRNA(seqData);
    // Convert data based for negative frame
    auto shift = frame;
    std::cerr << "Get ORF before " << l << " " << startLoc << " " << endLoc << std::endl;
    if (frame < 0)
    {
        shift = -frame;
        size_t org_end = endLoc;
        endLoc = l - startLoc - 1;
        startLoc = l - org_end - 1;
        //std::cout << startLoc << "||" << endLoc << std::endl;
        reverse(seqData);
        to35RNA(seqData);
    }
    std::cerr << "Get ORF after " << l << " " << startLoc << " " << endLoc << std::endl;
    shift -= 1;
    const std::set<std::string_view> endCodon{"UAA", "UAG", "UGA"};
    //                                         TTA    CTA    TCA
    const std::string startCodon = "AUG";
    //                              CAT
    std::vector<gene::GeneRange> result;
    std::string_view seqView(seqData.c_str(), l);
    
    #pragma omp parallel for shared(result, startCodon, endCodon, seqView)
    for (int64_t i = startLoc + shift; i < endLoc + shift; i += 3)
    {
        // Because of OpenMP, put the i+3 judge inside of loop
        if (i  < l - 3)
            // Check if current codon is start codon
            if (seqView.substr(i, 3).compare(startCodon) == 0)
                // Find if it has a end codon
                for (size_t j = i + 3; j < l - 3; j += 3)
                {
                    if (endCodon.find(seqView.substr(j, 3)) != endCodon.end())
                    {
                        auto start = (size_t) i;
                        auto end = (size_t) j + 2;
                        if (frame < 0)
                        {
                            start = l - start - 1;
                            end = l - end - 1;
                        }
                        #pragma omp critical
                        result.push_back({start, end, frame});
                        break;
                    }
                }
    }
    return result;
}