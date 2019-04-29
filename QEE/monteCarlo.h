#ifndef MONTECARLO_H
#define MONTECARLO_H


#include <string>

#include "tools.h"  // Numerical tools
#include "mtrnd.h"  // Random number generator


// Error measures
double epsilonP(const Vector &, const Vector &, int, double, const Vector &, double, double);
double epsilonP2(const Vector &, const Vector &, int, double, const Vector &, double);
double epsilonC(const Vector &, const Vector &, int, double, const Vector &, double);
double getMaxError(const Vector &, const Vector &, const IntVector &, double);
double getMaxError(const Vector &, const Vector &, double, const Vector &, double, double);

// MC functions
void setEnergies(const Vector &, const IntVector &, const std::vector<int> &, const std::vector<int> &, Vector &, std::vector<double> &);
void updateEnergies(const Vector &, const IntVector &, const std::vector<int> &, const std::vector<int> &, int, int, Vector &, std::vector<double> &);
void dynamicsRF(const Vector &, const IntVector &, int, double, std::vector<int> &, std::vector<int> &, Vector &, std::vector<double> &);
void updateCorrelations(const std::vector<int> &, const std::vector<int> &, Vector &);
void getSample(const Vector &, const IntVector &, MT::MersenneTwist &, int, int, int, std::vector<int> &, std::vector<int> &, Vector &, Vector &, std::vector<double> &);
void getFixedSample(const Vector &, const IntVector &, const std::vector<int> &, MT::MersenneTwist &, int, int, int, std::vector<int> &, std::vector<int> &, Vector &, Vector &, std::vector<double> &);
void writeFixedSample(FILE *, const Vector &, const IntVector &, const std::vector<int> &, MT::MersenneTwist &, int, int, int, std::vector<int> &, std::vector<int> &, Vector &, Vector &, std::vector<double> &);

// Auxiliary
void getNeff(const Vector &, int, int &, int &);
void getNeighbors(const Vector &, int, double, IntVector &);

// Control functions
void getError(const Vector &, const Vector &, int, double, int, int, double, double, std::vector<double> &);
void getErrorMCLearn(const Vector &, const Vector &, Vector &, double, int, int, double, double, Vector &, std::vector<double> &, std::vector<int> &);
void getEpitopeCorrelations(const Vector &, const Vector &, const std::vector<std::vector<int> > &, const std::vector<int> &, int, int, Vector &, std::vector<int> &);
void getEpitopeStates(FILE *, const Vector &, const Vector &, const std::vector<std::vector<int> > &, const std::vector<int> &, int, int, Vector &, std::vector<int> &);


#endif
