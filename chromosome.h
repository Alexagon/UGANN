#ifndef CHROMOSOME_H
#define CHROMOSOME_H

#include <QString>

struct chromosome {
    QString seq;
    double error;

    chromosome(QString s);

    chromosome(QString s, double e);

    bool operator <(const chromosome &x) const;
};

#endif // CHROMOSOME_H
