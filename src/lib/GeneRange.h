#pragma once
#ifndef _GENERANGE_H
#define _GENERANGE_H
#include <stdint.h>
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
}
#endif