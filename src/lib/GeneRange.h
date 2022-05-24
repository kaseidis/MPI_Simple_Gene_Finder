#pragma once
#ifndef _GENERANGE_H
#define _GENERANGE_H
#include <stdint.h>

#define INVALID_RANGE_LOC (~((unsigned long long) 0))
#define INVALID_FRAME (0)

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
         *        INVALID_RANGE_LOC indicates the object is invalid
         */
        unsigned long long start;
        /**
         * @brief end position of gene
         *        INVALID_RANGE_LOC indicates the object is invalid
         */
        unsigned long long end;
        /**
         * @brief frame position of gene
         *        INVALID_FRAME indicates the object is invalid
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

        /**
         * @brief Check if the range object is valid.
         *
         * @return true     It is a valid range object.
         * @return false    It is a invalid range object.
         */
        explicit operator bool() const {
            return (start!=INVALID_RANGE_LOC) && (end!=INVALID_RANGE_LOC) && (frame!=INVALID_FRAME);
        };
    } GeneRange;
}
#endif