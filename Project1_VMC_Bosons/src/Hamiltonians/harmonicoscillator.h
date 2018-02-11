#pragma once
#include "hamiltonian.h"
#include <vector>

class HarmonicOscillator : public Hamiltonian {
public:
    HarmonicOscillator(System* system, double omega,bool analytical);
    double computeLocalEnergy(std::vector<Particle*> particles);

    double computePotentialEnergy(std::vector<Particle*> particles);
    double computeKineticEnergy(std::vector<Particle*> particles);

private:
    double m_omega = 0;
    bool   m_analytical = true;
};

