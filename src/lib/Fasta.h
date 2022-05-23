#pragma once

#ifndef _FASTA_H
#define _FASTA_H

#include <iostream>
#include <string>
#include <fstream>

/**
 * @brief A fasta sequence object, contains a label and sequence.
 */
class Sequence
{
private:
    std::string label;
    std::string sequence;
    bool error;

public:
    /**
     * @brief Construct a new Sequence object
     *
     * @param label
     * @param seq
     */
    Sequence(const std::string &label, const std::string &seq);
    /**
     * @brief Construct a new Sequence object, which contains
     *        a validation flag.
     *
     * @param error True, if the object is not valid
     */
    Sequence(bool error = true);
    /**
     * @brief Get the Label object
     *
     * @return const std::string&
     */
    const std::string &getLabel() const;
    /**
     * @brief Get the Sequence object
     *
     * @return const std::string&
     */
    const std::string &getSequence() const;
    /**
     * @brief Set the Label object
     *
     * @param label
     */
    void setLabel(const std::string &label);
    /**
     * @brief Set the Sequence object
     *
     * @param seq
     */
    void setSequence(const std::string &seq);
    /**
     * @brief Check if the sequence object is valid.
     *
     * @return true     It is a valid sequence object.
     * @return false    It is a invalid sequence object.
     */
    explicit operator bool() const;
};

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