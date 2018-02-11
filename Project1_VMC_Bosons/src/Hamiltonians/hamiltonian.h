#pragma once
#include <vector>

class Hamiltonian {
public:
    Hamiltonian(class System* system);
    virtual double computeLocalEnergy(std::vector<class Particle*> particles) = 0;
    double computeNumericalKineticEnergy(std::vector<class Particle*> particles);
    double computeKineticEnergy(std::vector<class Particle*> particles);
    virtual double getKineticEnergy(std::vector<class Particle*> particles);
    virtual double getPotentialEnergy(std::vector<class Particle*> particles);
protected:
    class System* m_system = nullptr;
};

