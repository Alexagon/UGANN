#include "chromosome.h"

chromosome::chromosome(QString s)
{
    this->seq = s;
    this->error = 1;
}

chromosome::chromosome(QString s, double e)
{
    this->seq = s;
    this->error = e;
}

bool chromosome::operator <(const chromosome &x) const
{
    return (this->error < x.error);
}
