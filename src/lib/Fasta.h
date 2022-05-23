#pragma once

#ifndef _FASTA_H
#define _FASTA_H

#include <iostream>
#include <string>
#include <fstream>
#include "Sequence.h"

/**
 * @brief  A fasta file object, that parse fasta file.
 */
class Fasta
{
private:
    std::fstream file;
    std::string filename;
    std::string label;

public:
    /**
     * @brief Construct a new Fasta object
     *
     * @param filename
     * @param mode
     */
    Fasta(const char *filename, std::ios_base::openmode mode);
    /**
     * @brief close the file pointer of this Fasta object
     */
    void close();
    /**
     * @brief Destroy the Fasta object, do same operation
     *        as close()
     */
    ~Fasta();
    /**
     * @brief Write a sequence to file
     *
     * @param seq
     * @param lineWidth
     * @return true         Operation sucessful.
     * @return false        Operation failed.
     */
    bool write(const Sequence &seq, size_t lineWidth = 70);
    /**
     * @brief Get the Next Sequence object
     *
     * @return Sequence
     */
    Sequence getNextSequence();
};

#endif