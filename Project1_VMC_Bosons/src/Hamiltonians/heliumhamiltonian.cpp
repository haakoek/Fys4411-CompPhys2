#include "heliumhamiltonian.h"
#include "harmonicoscillator.h"
#include <cassert>
#include <iostream>
#include "../system.h"
#include "../particle.h"
#include "../WaveFunctions/wavefunction.h"
#include <cmath>

HeliumHamiltonian::HeliumHamiltonian(System* system, bool interaction) :
    Hamiltonian(system) {
        m_interaction = interaction;
}

HeliumHamiltonian::HeliumHamiltonian(System* system) :

    HeliumHamiltonian(system,true) {
}

double HeliumHamiltonian::computeLocalEnergy(std::vector<Particle *> particles) {

    double r12     = 0;
    double r1      = 0;
    double r2      = 0;
    double r1dotr2 = 0;

    double alpha = m_system->getWaveFunction()->getParameters()[0];
    double beta  = m_system->getWaveFunction()->getParameters()[1];
    double Z     = 2.0;

    for (int k=0; k<m_system->getNumberOfDimensions(); k++) {

        const double x1     = particles[0]->getPosition()[k];
        const double x2     = particles[1]->getPosition()[k];
        const double x12    = x2 - x1;
        r1  += x1  * x1;
        r2  += x2  * x2;
        r12 += x12 * x12;
        r1dotr2 += x1*x2;

    }

    r1 = sqrt(r1); r2 = sqrt(r2); r12 = sqrt(r12);

    const double localEnergy = (alpha-Z)*((1.0/r1) + (1.0/r2)) - alpha*alpha
                               + (1.0/(2.0*(1+beta*r12)*(1+beta*r12)))*((alpha/r12)*((r1+r2)*(1-r1dotr2/(r1*r2)))
                               - (2.0/r12) + (2.0*beta/(1+beta*r12)) - 1.0/(2.0*(1+beta*r12)*(1+beta*r12)));


    const double interactionEnergy  = (m_interaction ? 1.0 / r12 : 0.0);
    //const double kineticEnergy      = Hamiltonian::computeKineticEnergy(particles);
    //const double potentialEnergy    = -2.0 / sqrt(r1) - 2.0 / sqrt(r2);

    //return kineticEnergy + potentialEnergy + interactionEnergy;


    return localEnergy+interactionEnergy;
}

