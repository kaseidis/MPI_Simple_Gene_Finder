#include "Fasta.h"
#include <sstream>
#include <string>
#include <algorithm>
#include <omp.h>


Fasta::Fasta(const char *filename, std::ios_base::openmode mode)
{
    this->filename = filename;
    this->file = std::fstream(filename, mode);
    // Get first sequence label
    if ((mode & std::ios::in) != 0)
        while (std::getline(this->file, this->label))
        {
            // Trim CRLF
            if (this->label.length() != 0 &&
                this->label[this->label.length() - 1] == '\r')
            {
                this->label = this->label.substr(0, this->label.length() - 1);
            }
            // Find first label line
            if (this->label.length() != 0 &&
                this->label[0] == '>')
            {
                this->label = this->label.substr(1, this->label.length() - 1);
                break;
            }
        }
}

Fasta::~Fasta()
{
    this->close();
}

void Fasta::close()
{
    if (this->file.is_open())
        this->file.close();
}

bool Fasta::write(const Sequence &seq, size_t lineWidth)
{
    // Check file
    if (!this->file.is_open())
        return false;
    // Print label
    this->file << '>' << seq.getLabel() << std::endl;
    // Print sequence
    auto sequence = seq.getSequence();
    for (int64_t i = 0; i < sequence.length(); i += lineWidth)
        if (!(this->file << sequence.substr(i, lineWidth) << std::endl))
            return false;
    // Return sucessful
    return true;
}

Sequence Fasta::getNextSequence()
{
    // File not open
    if (!this->file.is_open())
    {
        return Sequence(true);
    }
    // File is eof
    if (this->file.eof())
    {
        return Sequence(true);
    }
    std::ostringstream seq;
    std::string line;
    while (std::getline(this->file, line))
    {
        // Trim CRLF
        if (line.length() != 0 &&
            line[line.length() - 1] == '\r')
        {
            line = line.substr(0, line.length() - 1);
        }
        // If got new label line, return last sequence
        if (line[0] == '>')
        {
            auto result = Sequence(this->label, seq.str());
            this->label = line.substr(1, line.length() - 1);
            return result;
        }
        // Standarize sequence
        #pragma omp parallel for
        for (int i = 0; i < line.length(); ++i)
        {
            line[i] = line[i] == '_' ? '-' : (char)std::toupper(line[i]);
        }
        // Build sequence
        seq << line;
    }
    // Deal with last sequence
    if (seq.str().length() != 0 || this->label.length() != 0)
    {
        auto result = Sequence(this->label, seq.str());
        this->label = "";
        return result;
    }
    // File is eof
    return Sequence(true);
}