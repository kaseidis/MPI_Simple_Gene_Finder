#include "./lib/Fasta.h"
#include "./lib/InputParser.h"
#include "./lib/orf_finder.h"
#include "./lib/gene_judge.h"
#include <iostream>
#include <vector>
#include <omp.h>
#include <string>
#include <atomic>
#include <memory>
#include <sstream>
#include <chrono>

/**
 * @brief C++11 version of sprintf
 *        Reference: https://stackoverflow.com/questions/2342162/
 * @tparam Args
 * @param format
 * @param args
 * @return std::string
 */
template <typename... Args>
std::string string_format(const std::string &format, Args... args)
{
    int size_s = std::snprintf(nullptr, 0, format.c_str(), args...) + 1; // Extra space for '\0'
    if (size_s <= 0)
    {
        throw std::runtime_error("Error during formatting.");
    }
    auto size = static_cast<size_t>(size_s);
    std::unique_ptr<char[]> buf(new char[size]);
    std::snprintf(buf.get(), size, format.c_str(), args...);
    return std::string(buf.get(), buf.get() + size - 1); // We don't want the '\0' inside
}

/**
 * @brief Get the gene object, take isGene function from dynamic linked lib.
 *
 * @param orfs    vector that contains ORFS, to check if it is a gene.
 * @param seq     Sequence to judge.
 * @param start   Index of start ORF object
 * @param end     Index of end ORF object.
 * @return std::vector<gene::GeneRange>
 */
std::vector<gene::GeneRange> get_gene(
    const std::vector<gene::GeneRange> &orfs, const Sequence &seq,
    size_t start, size_t end)
{
    std::vector<gene::GeneRange> result;
    result.resize(orfs.size());
    std::atomic<size_t> resultIndex = 0;
    #pragma omp parallel for
    for (int64_t i = start; i < end; ++i)
    {
        auto resSeq = isGene(orfs[i], seq);
        if (resSeq)
            result[resultIndex++] = resSeq;
    }
    result.resize(resultIndex);
    return result;
}

/**
 * @brief Finding gene from fasta and save it to another fasta file.
 * 
 * @param input_filepath 
 * @param output_filepath 
 * @param print_pattern 
 * @param line_width 
 * @return int 
 */
int finding_gene(const char *input_filepath, const char *output_filepath,
         const char *print_pattern, size_t line_width = 70)
{
    // Open files
    Fasta f(input_filepath, std::ios::in);
    Fasta f_out(output_filepath, std::ios::out);
    // Get all sequences
    for (auto seq = f.getNextSequence(); seq; seq = f.getNextSequence())
    {
        
        for (int frame = -3; frame <= 3; ++frame) {
            if (frame==0)
                continue;
            // Get orfs
            auto orfs = gene::getORFS(seq, frame, 0,
                                    seq.getSequence().length());
            // Filter orfs
            auto g = get_gene(orfs, seq, 0, orfs.size());
            // Save gene to file
            for (auto i = 0; i < g.size(); i++)
            {
                Sequence seq_out(
                    string_format(print_pattern,seq.getLabel().c_str(),frame,g[i].start,g[i].end),
                    seq.getSequence().substr(
                        g[i].abs_start(), g[i].length()));
                f_out.write(seq_out, line_width);
            }
        }
    }
    // Close file
    f.close();
    f_out.close();
    // Return 0 for sucessful.
    return 0;
}

/**
 * @brief Print usage of program
 * 
 * @param prog program name
 */
void print_usage(const char* prog)
{
    std::cout << "Usage: " << prog << " --input INPUT_FILE_PATH"
              << " --output OUTPUT_FILE_PATH"
              << " [--pattern LABEL_PATTERN --output-line-width WIDTH --time]" << std::endl;
    std::cout << "    Default:" << std::endl <<
        "        LABEL_PATTERN = '%s | gene | frame=%d | LOC=[%d,%d]'" << std::endl <<
        "        WIDTH = 70" << std::endl;
}

int main(int argc, char **argv)
{
    InputParser input = InputParser(argc, argv);

    // Check for helper option
    if (input.cmdOptionExists("-h") || input.cmdOptionExists("--help")) {
        print_usage(argv[0]);
        return 1;
    }

    // Check for input and output option
    std::string input_file, output_file, pattern = "%s | gene | frame=%d | LOC=[%d,%d]";
    size_t line_width = 70;
    if (input.cmdOptionExists("--input") && input.cmdOptionExists("--output"))
    {
        input_file = input.getCmdOption("--input");
        output_file = input.getCmdOption("--output");
    }
    else
    {
        std::cerr << "Invalid argument" << std::endl;
        print_usage(argv[0]);
        return 1;
    }
    // Check for pattern option
    if (input.cmdOptionExists("--pattern"))
        pattern = input.getCmdOption("--pattern");
    // check for --output-line-width option
    if (input.cmdOptionExists("--output-line-width"))
    {
        auto line_width_option = input.getCmdOption("--output-line-width");
        std::istringstream line_width_stream(line_width_option);
        line_width_stream >> line_width;
    }
    // check for --time option
    bool check_time = false;
    if (input.cmdOptionExists("--time"))
    {
        auto line_width_option = input.getCmdOption("--time");
        check_time = true;
    }
    auto start = std::chrono::high_resolution_clock::now();
    auto result = finding_gene(input_file.c_str(), output_file.c_str(), pattern.c_str(),line_width);
    // Timing
    if (check_time) {
        auto finish = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = finish - start;
        std::cout << elapsed.count() << std::endl;
    }
    return result;
}