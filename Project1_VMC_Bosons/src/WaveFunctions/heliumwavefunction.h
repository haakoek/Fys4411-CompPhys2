#pragma once
#include "wavefunction.h"

class HeliumWaveFunction : public WaveFunction
{
public:
    HeliumWaveFunction(class System* system, double alpha, double beta);
    double evaluate(std::vector<Particle *> particles);
    std::vector<double> computeGradient(std::vector<class Particle *> particles);
    double computeAlphaDerivative(std::vector<Particle *> particles);
};


