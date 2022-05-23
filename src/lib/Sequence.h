#pragma once
#ifndef _SEQUENCE_H
#define _SEQUENCE_H
#include <string>
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
#endif