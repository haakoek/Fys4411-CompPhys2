#include "hamiltonian.h"
#include "../system.h"
#include "../particle.h"
#include "../WaveFunctions/wavefunction.h"
#include "../WaveFunctions/simplegaussian.h"

Hamiltonian::Hamiltonian(System* system) {
    m_system = system;
}


double Hamiltonian::computeKineticEnergy(std::vector<class Particle*> particles) {
    return computeNumericalKineticEnergy(particles);
}

double Hamiltonian::computeNumericalKineticEnergy(std::vector<Particle*> particles) {

    double h = 1e-6;
    double laplacian = 0;

    double psi = m_system->getWaveFunction()->evaluate(particles);
    double h2  = h*h;

    for(int p = 0; p < m_system->getNumberOfParticles(); p++) {

            for(int d = 0; d < m_system->getNumberOfDimensions(); d++) {

                particles[p]->adjustPosition(h,d);

                double psi_p = m_system->getWaveFunction()->evaluate(particles);

                particles[p]->adjustPosition(-2.0*h,d);
                double psi_m = m_system->getWaveFunction()->evaluate(particles);

                particles[p]->adjustPosition(h,d);


                laplacian += (psi_p - 2.0*psi + psi_m)/(h2);

            }

    }

    return -0.5*laplacian/psi;

}

double Hamiltonian::getKineticEnergy(std::vector<Particle *> particles) {
    return 0.0;
}

double Hamiltonian::getPotentialEnergy(std::vector<Particle *> particles) {
    return 0.0;
}
