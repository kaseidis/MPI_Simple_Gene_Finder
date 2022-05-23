#include "./lib/Fasta.h"
#include "./lib/InputParser.h"
#include "./lib/orf_finder.h"
#include <iostream>
#include <vector>
#include "gene_judge.h"
#include <omp.h>
#include <string>
#include <atomic>
#include <memory>
#include <sstream>
#include <mpi.h>

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
    const std::vector<gene::GeneRange> &orfs,
    const Sequence &seq, size_t start, size_t end)
{
    std::vector<gene::GeneRange> result;
    result.resize(orfs.size());
    std::atomic<size_t> resultIndex = 0;
#pragma omp parallel for
    for (int64_t i = start; i < end; ++i)
    {
        if (isGene(orfs[i], seq))
            result[resultIndex++] = orfs[i];
    }
    result.resize(resultIndex);
    return result;
}

inline size_t get_job_count(size_t total_job_count, size_t mpi_rank, size_t mpi_size)
{
    auto additional = mpi_rank < total_job_count % mpi_size ? 1 : 0;
    return total_job_count / mpi_size + additional;
}

int get_need(unsigned long long *job_counts, size_t total_job, size_t mpi_size)
{
    for (auto i = 0; i < mpi_size; ++i)
        if (job_counts[i] < get_job_count(total_job, i, mpi_size))
            return i;
    return -1;
}

int get_full(unsigned long long *job_counts, size_t total_job, size_t mpi_size)
{
    for (auto i = 0; i < mpi_size; ++i)
        if (job_counts[i] > get_job_count(total_job, i, mpi_size))
            return i;
    return -1;
}

void send_gene_range(std::vector<gene::GeneRange> &ranges, MPI_Datatype MPI_GENE_TYPE, size_t count, int target)
{
    gene::GeneRange range_buf[count];
    for (size_t i = 0; i < count; ++i)
    {
        range_buf[i] = ranges.back();
        ranges.pop_back();
    }
    MPI_Send(range_buf, count, MPI_GENE_TYPE, target, 0, MPI_COMM_WORLD);
}

MPI_Status recv_gene_range(std::vector<gene::GeneRange> &ranges, MPI_Datatype MPI_GENE_TYPE, size_t count, int target)
{
    gene::GeneRange range_buf[count];
    MPI_Status s;
    MPI_Recv(range_buf, count, MPI_GENE_TYPE, 0, 0, MPI_COMM_WORLD, &s);
    for (size_t i = 0; i < count; ++i)
        ranges.push_back(range_buf[i]);
    return s;
}

int findingGene(const char *input_filepath, const char *output_filepath,
                const char *print_pattern, int mpi_rank, int mpi_size, size_t line_width = 70)
{
    // Create type for gene range
    const int nitems = 3;
    int blocklengths[3] = {1, 1, 1};
    MPI_Datatype types[3] = {MPI_UNSIGNED_LONG_LONG, MPI_UNSIGNED_LONG_LONG, MPI_UINT8_T};
    MPI_Datatype MPI_GENE_RANGE;
    MPI_Aint offsets[3];
    offsets[0] = offsetof(gene::GeneRange, start);
    offsets[1] = offsetof(gene::GeneRange, end);
    offsets[2] = offsetof(gene::GeneRange, frame);
    MPI_Type_create_struct(nitems, blocklengths, offsets, types, &MPI_GENE_RANGE);
    MPI_Type_commit(&MPI_GENE_RANGE);

    // Reading orfs from file
    std::vector<gene::GeneRange> local_orfs;
    std::vector<gene::GeneRange> local_gene;
    std::vector<gene::GeneRange> all_gene;
    Fasta f(input_filepath, std::ios::in);
    for (auto seq = f.getNextSequence(); seq; seq = f.getNextSequence())
    {
        auto orfs = gene::getORFS(seq, -2, 0,
                                  seq.getSequence().length());
        // Store result to local orfs vector
        if (local_orfs.capacity() < local_orfs.size() + orfs.size())
        {
            local_orfs.reserve(local_orfs.size() + orfs.size());
        }
        local_orfs.insert(local_orfs.end(), orfs.begin(), orfs.end());
        // Balancing ORFS
        unsigned long long job_count = local_orfs.size();
        MPI_Allreduce(&job_count, &job_count, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
        // If is main process allocate job
        if (mpi_rank == 0)
        {
            // Get job count from every sub node
            unsigned long long job_counts[mpi_size];
            job_counts[0] = local_orfs.size();
            for (auto i = 1; i < mpi_size; ++i)
            {
                MPI_Status s;
                MPI_Recv(&(job_counts[i]), 1, MPI_UNSIGNED_LONG_LONG, i, 0, MPI_COMM_WORLD, &s);
            }
            // Balancing nodes
            auto send_node = get_full(job_counts, job_count, mpi_size);
            auto recv_node = get_need(job_counts, job_count, mpi_size);
            while (send_node != -1 || recv_node != -1)
            {
                send_node = get_full(job_counts, job_count, mpi_size);
                recv_node = get_need(job_counts, job_count, mpi_size);
                size_t send_count = job_counts[send_node] - get_job_count(job_count, send_node, mpi_size);
                size_t recv_count = get_job_count(job_count, recv_node, mpi_size) - job_counts[recv_node];
                unsigned long long current_count = send_count < recv_count ? send_count : recv_count;
                // Send count first and then send target node
                if (send_node != 0)
                {
                    MPI_Send(&current_count, 1, MPI_UNSIGNED_LONG_LONG, send_node, 0, MPI_COMM_WORLD);
                    MPI_Send(&recv_node, 1, MPI_INT, send_node, 0, MPI_COMM_WORLD);
                }
                if (recv_node != 0)
                {
                    MPI_Send(&current_count, 1, MPI_UNSIGNED_LONG_LONG, recv_node, 0, MPI_COMM_WORLD);
                    MPI_Send(&send_node, 1, MPI_INT, recv_node, 0, MPI_COMM_WORLD);
                }
                // Record balancing result
                job_counts[send_node] -= current_count;
                job_counts[recv_node] += current_count;

                // If main node need to do job, do it right now
                if (send_node == 0)
                    send_gene_range(local_orfs, MPI_GENE_RANGE, current_count, recv_node);
                if (recv_node == 0)
                    recv_gene_range(local_orfs, MPI_GENE_RANGE, current_count, send_node);
            }
        }
        else // If it is sub node
        {
            // Send current job count to main node for balancing
            unsigned long long local_size = local_orfs.size();
            MPI_Send(&local_size, 1, MPI_UNSIGNED_LONG_LONG, 0, 0, MPI_COMM_WORLD);

            // Calculate target job count
            auto addition_rank_max = job_count % mpi_size;
            job_count = get_job_count(job_count, mpi_rank, mpi_size);

            // Getting balancing target from main node. Then, getting job from sub node.
            while (local_orfs.size() != job_count)
            {
                unsigned long long current_count;
                int target_node;
                MPI_Status s;
                MPI_Recv(&current_count, 1, MPI_UNSIGNED_LONG_LONG, 0, 0, MPI_COMM_WORLD, &s);
                MPI_Recv(&target_node, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &s);
                if (local_orfs.size() > job_count)
                    send_gene_range(local_orfs, MPI_GENE_RANGE, current_count, target_node);
                else
                    recv_gene_range(local_orfs, MPI_GENE_RANGE, current_count, target_node);
            }
        }

        // Getting gene
        auto gene_result = get_gene(local_orfs, seq, 0, orfs.size());

        // Get total gene count
        job_count = gene_result.size();
        MPI_Reduce(&job_count, &job_count, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);

        // Gethering gene to main node
        gene::GeneRange *final_result_buffer = NULL;
        if (mpi_rank == 0)
            final_result_buffer = (gene::GeneRange *)malloc(sizeof(gene::GeneRange) * job_count);
        MPI_Gather(gene_result.data(), gene_result.size(),
                   MPI_GENE_RANGE, final_result_buffer, job_count,
                   MPI_GENE_RANGE, 0, MPI_COMM_WORLD);
        // If it is main node, save result to file
        if (mpi_rank == 0)
        {
            Fasta f_out(output_filepath, std::ios::out | std::ios::app);
            for (unsigned long long i = 0; i < job_count; i++)
            {
                Sequence seq_out(
                    string_format(print_pattern, seq.getLabel().c_str(),
                                  final_result_buffer[i].start,
                                  final_result_buffer[i].end),
                    seq.getSequence().substr(
                        final_result_buffer[i].abs_start(),
                        final_result_buffer[i].length()));
                f_out.write(seq_out, line_width);
            }
            f_out.close();
        }
    }
    f.close();

    return 0;
}

void print_usage(const char *prog)
{
    std::cout << "Usage: " << prog << " --input INPUT_FILE_PATH"
              << " --output OUTPUT_FILE_PATH"
              << " [--pattern LABEL_PATTERN --output-line-width WIDTH]" << std::endl;
    std::cout << "    Default:" << std::endl
              << "        LABEL_PATTERN = '%s | gene | LOC=[%d,%d]'" << std::endl
              << "        WIDTH = 70" << std::endl;
}

int main(int argc, char **argv)
{

    // Init MPI
    int rank, size;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    InputParser input = InputParser(argc, argv);

    // Check for helper option
    if (input.cmdOptionExists("-h") || input.cmdOptionExists("--help"))
    {
        print_usage(argv[0]);
        return 1;
    }

    // Check for input and output option
    std::string input_file, output_file, pattern = "%s | gene | LOC=[%d,%d]";
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
    auto result = findingGene(input_file.c_str(), output_file.c_str(), pattern.c_str(), line_width, rank, size);
    MPI_Finalize();
    return result;
}