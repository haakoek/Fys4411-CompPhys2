#include "wavefunction.h"
#include "../particle.h"
#include <vector>
#include "../system.h"
#include <iostream>
using std::cout;
using std::endl;

WaveFunction::WaveFunction(System* system) {
    m_system = system;
}

std::vector<double> WaveFunction::computeGradient(std::vector<Particle*> particles) {

    cout << "balle" << endl;

    std::vector<double> gradient(m_system->getNumberOfDimensions()*m_system->getNumberOfParticles());
    double h = 1e-5;

    for(int p = 0; p < m_system->getNumberOfParticles(); p++) {

        for(int d = 0; d < m_system->getNumberOfDimensions(); d++) {

            particles[p]->adjustPosition(h,d);
            double psi_p = evaluate(particles);
            particles[p]->adjustPosition(-2*h,d);
            double psi_m = evaluate(particles);
            particles[p]->adjustPosition(h,d);

            gradient[m_system->getNumberOfDimensions()*p+d] = (psi_p-psi_m)/(2.0*h);


        }

    }

    return gradient;
}

void WaveFunction::setAlpha(double alpha) {
    m_parameters[1] = alpha;
}

void WaveFunction::setBeta(double beta) {
    m_parameters[2] = beta;
}

double WaveFunction::computeBetaDerivative(std::vector<Particle *> particles) {
    return 1;
}
