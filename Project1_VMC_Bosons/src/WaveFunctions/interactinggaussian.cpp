#include "interactinggaussian.h"
#include "simplegaussian.h"
#include <cmath>
#include <cassert>
#include "wavefunction.h"
#include "../system.h"
#include "../particle.h"

#include <iostream>

using namespace std;

InteractingGaussian::InteractingGaussian(System* system, double alpha, double beta) :
    WaveFunction(system) {
    assert(alpha >= 0);
    m_numberOfParameters = 2;
    m_parameters.reserve(2);
    m_parameters.push_back(alpha);
    m_parameters.push_back(beta);

}

double InteractingGaussian::evaluate(std::vector<Particle *> particles) {

    double psi = 1.0;
    const double alpha = m_parameters[0];
    const double beta  = m_parameters[1];

    for(int i = 0; i < m_system->getNumberOfParticles(); i++) {

        double r2 = 0;
        for(int j=0; j < m_system->getNumberOfDimensions(); j++) {

              if(j < 2) {
                r2+= particles[i]->getPosition()[j]*particles[i]->getPosition()[j];
              } else {
                r2+= beta*particles[i]->getPosition()[j]*particles[i]->getPosition()[j];
              }

        }
        psi*=exp(-alpha*r2);
    }

    psi*=correlation(particles); // regn ut prod_ij f(i,j)

    return psi;

}

double InteractingGaussian::correlation(std::vector<Particle *> particles) {

    double corr = 1.0;

    for(int i = 0; i < m_system->getNumberOfParticles(); i++) {
        for(int j = i+1; j < m_system->getNumberOfParticles(); j++) {
            corr*=f(particles[i],particles[j]);
        }
    }

    return corr;
}

double InteractingGaussian::f(Particle* pi, Particle* pj) {

    const double a = 0.0043; //??? a/a_ho = 0.0043

    double rij = 0.0;

    for(int d = 0; d < m_system->getNumberOfDimensions(); d++) {
        rij+= (pi->getPosition()[d]-pj->getPosition()[d])*(pi->getPosition()[d]-pj->getPosition()[d]);
    }

    rij = sqrt(rij);

    if(rij <= a) {
        return 0.0;
    } else {
        return 1.0-(a/rij);
    }

}

std::vector<double> InteractingGaussian::computeGradient(std::vector<Particle *> particles) {

    std::vector<double> gradient(m_system->getNumberOfDimensions()*m_system->getNumberOfParticles());
    double psi   = evaluate(particles);
    double alpha = m_parameters[0];
    double beta  = m_parameters[1];
    double a = 0.0043;

    for(int k = 0; k < m_system->getNumberOfParticles(); k++) {

        for(int d = 0; d < m_system->getNumberOfDimensions(); d++) {

            double xk  = particles[k]->getPosition()[d];
            double sum = 0.0;

            for(int j = 0; j < m_system->getNumberOfParticles(); j++) {

                if(j != k) {
                    double xj = particles[j]->getPosition()[d];
                    double r_kj = 0.0;

                    for(int d2 = 0; d2 < m_system->getNumberOfDimensions(); d2++) {
                        r_kj+=(particles[k]->getPosition()[d2] - particles[j]->getPosition()[d2])*(particles[k]->getPosition()[d2] - particles[j]->getPosition()[d2]);
                    }



                    sum+=(xk-xj)/pow(r_kj,1.5);
                }
            }

            if(d < 2) {
                gradient[m_system->getNumberOfDimensions()*k+d] = (-2.0*alpha*xk + a*sum)*psi;
            } else {
                gradient[m_system->getNumberOfDimensions()*k+d] = (-2.0*alpha*beta*xk + a*sum)*psi;
            }

        }


    }


    return gradient;
}

double InteractingGaussian::computeAlphaDerivative(std::vector<Particle *> particles) {

    double r_sum = 0.0;
    double beta  = m_parameters[1];

    for(int p = 0; p < m_system->getNumberOfParticles(); p++) {
        for(int d = 0; d < m_system->getNumberOfDimensions(); d++) {

            if(d < 2) {

                r_sum+= particles[p]->getPosition()[d]*particles[p]->getPosition()[d];
            } else {
                r_sum+= beta*particles[p]->getPosition()[d]*particles[p]->getPosition()[d];
            }
        }
    }

    return -r_sum;//*evaluate(particles);
}
