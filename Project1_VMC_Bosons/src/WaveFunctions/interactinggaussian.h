#pragma once
#include "wavefunction.h"

class InteractingGaussian : public WaveFunction {
public:
    InteractingGaussian(class System* system, double alpha, double beta);
    double evaluate(std::vector<class Particle*> particles);
    double correlation(std::vector<class Particle *> particles);
    double f(class Particle* pi, class Particle* pj);
    double computeAlphaDerivative(std::vector<class Particle *> particles);
    std::vector<double> computeGradient(std::vector<Particle*> particles);
};


