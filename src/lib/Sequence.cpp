#include "Sequence.h"

Sequence::Sequence(const std::string &label, const std::string &seq)
{
    this->label = label;
    this->sequence = seq;
    this->error = false;
}

Sequence::Sequence(bool error)
{
    this->error = error;
}

const std::string &Sequence::getLabel() const
{
    return this->label;
}

const std::string &Sequence::getSequence() const
{
    return this->sequence;
}

void Sequence::setLabel(const std::string &label)
{
    this->label = label;
}

void Sequence::setSequence(const std::string &seq)
{
    this->sequence = seq;
}

Sequence::operator bool() const
{
    return !this->error;
}