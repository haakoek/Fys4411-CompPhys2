#pragma once
#include "hamiltonian.h"
#include <vector>

class HeliumHamiltonian : public Hamiltonian
{
public:
    HeliumHamiltonian(class System* system);
    HeliumHamiltonian(class System* system, bool interaction);
    double computeLocalEnergy(std::vector<Particle *> particles);
private:
    bool m_interaction = false;
};
