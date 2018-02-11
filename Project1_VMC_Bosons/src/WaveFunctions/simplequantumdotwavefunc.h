#pragma once
#include "wavefunction.h"

class SimpleQuantumDotWaveFunc : public WaveFunction
{
public:
    SimpleQuantumDotWaveFunc(class System* system, double omega, double alpha, double beta);
    double evaluate(std::vector<Particle *> particles);
    double computeAlphaDerivative(std::vector<Particle *> particles);
    double computeBetaDerivative(std::vector<Particle *> particles);
    std::vector<double> computeGradient(std::vector<Particle *> particles);
};

