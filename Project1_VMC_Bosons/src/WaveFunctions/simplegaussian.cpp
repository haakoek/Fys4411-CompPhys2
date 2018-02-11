#include "simplegaussian.h"
#include <cmath>
#include <cassert>
#include "wavefunction.h"
#include "../system.h"
#include "../particle.h"

#include <iostream>

using std::cout;
using std::endl;

SimpleGaussian::SimpleGaussian(System* system, double alpha) :
        WaveFunction(system) {
    assert(alpha >= 0);
    m_numberOfParameters = 1;
    m_parameters.reserve(1);
    m_parameters.push_back(alpha);
}

double SimpleGaussian::evaluate(std::vector<class Particle*> particles) {
    /* You need to implement a Gaussian wave function here. The positions of
     * the particles are accessible through the particle[i].getPosition()
     * function.
     *
     * For the actual expression, use exp(-alpha * r^2), with alpha being the
     * (only) variational parameter.
     */

    double psi = 1.0;
    double alpha = m_parameters[0];

    for(int i = 0; i < m_system->getNumberOfParticles(); i++) {

        std::vector<double> ri = particles[i]->getPosition();
        double r2 = 0;

        for(int j=0; j < m_system->getNumberOfDimensions(); j++) {
              r2+= ri[j]*ri[j];
        }


        psi*=exp(-alpha*r2);
    }



    return psi;
}

std::vector<double> SimpleGaussian::computeGradient(std::vector<Particle *> particles) {
    std::vector<double> gradient(m_system->getNumberOfDimensions()*m_system->getNumberOfParticles());
    double alpha = m_parameters[0];

    for(int p = 0; p < m_system->getNumberOfParticles(); p++) {
        for(int d = 0; d < m_system->getNumberOfDimensions(); d++) {
            gradient[m_system->getNumberOfDimensions()*p+d] = -2.0*alpha*particles[p]->getPosition()[d]*evaluate(particles);
        }
    }

    return gradient;

}



double SimpleGaussian::computeAlphaDerivative(std::vector<Particle *> particles) {
    double r_sum = 0.0;
    for(int p = 0; p < m_system->getNumberOfParticles(); p++) {
        for(int d = 0; d < m_system->getNumberOfDimensions(); d++) {
            r_sum+= particles[p]->getPosition()[d]*particles[p]->getPosition()[d];
        }
    }

    return -r_sum;//*evaluate(particles);
}
