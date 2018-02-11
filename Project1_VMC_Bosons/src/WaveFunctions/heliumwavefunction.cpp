#include "heliumwavefunction.h"
#include "interactinggaussian.h"
#include "simplegaussian.h"
#include <cmath>
#include <cassert>
#include "wavefunction.h"
#include "../system.h"
#include "../particle.h"

#include <iostream>

using namespace std;

HeliumWaveFunction::HeliumWaveFunction(System* system, double alpha, double beta) :
    WaveFunction(system) {
        assert(alpha >= 0);
        assert(beta  >= 0);
        m_numberOfParameters = 2;
        m_parameters.reserve(2);
        m_parameters.push_back(alpha);
        m_parameters.push_back(beta);

}

double HeliumWaveFunction::evaluate(std::vector<Particle *> particles) {

    double a   = 0.5;
    double r1  = 0.0;
    double r2  = 0.0;
    double r12 = 0.0;

    double alpha = m_parameters[0]; double beta = m_parameters[1];

    for(int d = 0; d < m_system->getNumberOfDimensions(); d++) {
        double x1= particles[0]->getPosition()[d];
        double x2= particles[1]->getPosition()[d];
        r1+=x1*x1;
        r2+=x2*x2;
        r12+=(x1-x2)*(x1-x2);
    }

    r1  = sqrt(r1);
    r2  = sqrt(r2);
    r12 = sqrt(r12);

    return exp(-alpha*(r1+r2))*exp((a*r12)/(1.0+beta*r12));
}

std::vector<double> HeliumWaveFunction::computeGradient(std::vector<Particle *> particles) {


    double alpha = m_parameters[0];
    double beta  = m_parameters[1];
    std::vector<double> gradient(6);

    double a   = 0.5;
    double r1  = 0.0;
    double r2  = 0.0;
    double r12 = 0.0;

    for(int d = 0; d < m_system->getNumberOfDimensions(); d++) {
        double x1= particles[0]->getPosition()[d];
        double x2= particles[1]->getPosition()[d];
        r1+=x1*x1;
        r2+=x2*x2;
        r12+=(x1-x2)*(x1-x2);
    }

    r1  = sqrt(r1);
    r2  = sqrt(r2);
    r12 = sqrt(r12);

    double psi = evaluate(particles);

    double x1 = particles[0]->getPosition()[0];
    double y1 = particles[0]->getPosition()[1];
    double z1 = particles[0]->getPosition()[2];

    double x2 = particles[1]->getPosition()[0];
    double y2 = particles[1]->getPosition()[1];
    double z2 = particles[1]->getPosition()[2];

    gradient[0] = ((-alpha*x1/r1) + (a*(x1-x2)/r12)/((1+beta*r12)*(1+beta*r12)))*psi;
    gradient[1] = ((-alpha*y1/r1) + (a*(y1-y2)/r12)/((1+beta*r12)*(1+beta*r12)))*psi;
    gradient[2] = ((-alpha*z1/r1) + (a*(z1-z2)/r12)/((1+beta*r12)*(1+beta*r12)))*psi;

    gradient[3] = ((-alpha*x2/r2) - (a*(x1-x2)/r12)/((1+beta*r12)*(1+beta*r12)))*psi;
    gradient[4] = ((-alpha*y2/r2) - (a*(y1-y2)/r12)/((1+beta*r12)*(1+beta*r12)))*psi;
    gradient[5] = ((-alpha*z2/r2) + (a*(z1-z2)/r12)/((1+beta*r12)*(1+beta*r12)))*psi;


    return gradient;

}

double HeliumWaveFunction::computeAlphaDerivative(std::vector<Particle *> particles) {

    double r_sum = 0.0;

    for(int p = 0; p < m_system->getNumberOfParticles(); p++) {
        for(int d = 0; d < m_system->getNumberOfDimensions(); d++) {
                r_sum+= particles[p]->getPosition()[d]*particles[p]->getPosition()[d];
        }
    }

    return -r_sum*evaluate(particles);
}


