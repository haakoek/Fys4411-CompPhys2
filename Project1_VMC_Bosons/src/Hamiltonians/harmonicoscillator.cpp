#include "harmonicoscillator.h"
#include <cassert>
#include <iostream>
#include "../system.h"
#include "../particle.h"
#include "../WaveFunctions/wavefunction.h"

using std::cout;
using std::endl;

HarmonicOscillator::HarmonicOscillator(System* system, double omega, bool analytical) :
        Hamiltonian(system) {
    assert(omega > 0);
    m_omega  = omega;
    m_analytical = analytical;
}

double HarmonicOscillator::computeLocalEnergy(std::vector<Particle*> particles) {
    /* Here, you need to compute the kinetic and potential energies. Note that
     * when using numerical differentiation, the computation of the kinetic
     * energy becomes the same for all Hamiltonians, and thus the code for
     * doing this should be moved up to the super-class, Hamiltonian.
     *
     * You may access the wave function currently used through the
     * getWaveFunction method in the m_system object in the super-class, i.e.
     * m_system->getWaveFunction()...
     */

    double potentialEnergy   = computePotentialEnergy(particles);

    double kineticEnergy     = 0.0;


    if (m_analytical) {
        kineticEnergy = computeKineticEnergy(particles);
    } else {
        kineticEnergy = computeNumericalKineticEnergy(particles);
    }

    return kineticEnergy + potentialEnergy;
}


double HarmonicOscillator::computePotentialEnergy(std::vector<Particle*> particles) {

    double potentialEnergy = 0.0;

    for(int i = 0; i < m_system->getNumberOfParticles(); i++) {

        double r2 = 0;

        for(int j=0; j < m_system->getNumberOfDimensions(); j++) {
                r2 += particles[i]->getPosition()[j]*particles[i]->getPosition()[j];
        }
        potentialEnergy+=r2;
    }

    potentialEnergy = 0.5*m_omega*m_omega*potentialEnergy;

    return potentialEnergy;
}

double HarmonicOscillator::computeKineticEnergy(std::vector<Particle*> particles) {

    double alpha = m_system->getWaveFunction()->getParameters()[0];
    double d = m_system->getNumberOfDimensions();
    double N = m_system->getNumberOfParticles();

    double kineticEnergy = 0.0;

    for(int i = 0; i < m_system->getNumberOfParticles(); i++) {

        double r2 = 0;

        for(int j=0; j < m_system->getNumberOfDimensions(); j++) {
                r2 += particles[i]->getPosition()[j]*particles[i]->getPosition()[j];
        }

        kineticEnergy+=r2;
    }

    kineticEnergy = d*N*alpha - 2.0*alpha*alpha*kineticEnergy;

    return kineticEnergy;

}
