#pragma once
#include "wavefunction.h"

class SimpleGaussian : public WaveFunction {
public:
    SimpleGaussian(class System* system, double alpha);
    double evaluate(std::vector<class Particle*> particles);
    double computeAlphaDerivative(std::vector<Particle *> particles);
    std::vector<double> computeGradient(std::vector<Particle *> particles);
};
