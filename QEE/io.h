#ifndef IO_H
#define IO_H


#include <vector>
#include <math.h>
#include <stdio.h>

#include "tools.h"  // Numerical tools


void getCorrelations(FILE *, Vector &);
void getCouplings(FILE *, Vector &);
void getSS(FILE *, IntVector &p);
void printCouplings(FILE *, const Vector &);
void printSupplementaryOutput(FILE *, double, const std::vector<double> &, double, int, unsigned long, unsigned long);


#endif
