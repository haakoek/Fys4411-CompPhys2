#pragma once
#include "hamiltonian.h"

class EllipticInteractingOscillator : public Hamiltonian
{
public:
    EllipticInteractingOscillator(System* system, double gamma);
    double computeLocalEnergy(std::vector<class Particle *> particles);
    double computePotentialEnergy(std::vector<class Particle *> particles);
    double computeKineticEnergy(std::vector<class Particle *> particles);
    double interaction(Particle* pi, Particle* pj);


private:
    double m_gamma = 0;
};

