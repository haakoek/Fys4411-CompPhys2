#include "simplequantumdotwavefunc.h"
#include <cmath>
#include <cassert>
#include "wavefunction.h"
#include "../system.h"
#include "../particle.h"
#include <iomanip>
#include <iostream>

using namespace std;

SimpleQuantumDotWaveFunc::SimpleQuantumDotWaveFunc(System* system,double omega, double alpha, double beta) :
    WaveFunction(system) {
    assert(alpha >= 0);
    assert(beta  >= 0);
    assert(omega >= 0);
    m_numberOfParameters = 3;
    m_parameters.reserve(3);
    m_parameters.push_back(omega);
    m_parameters.push_back(alpha);
    m_parameters.push_back(beta);
}

double SimpleQuantumDotWaveFunc::evaluate(std::vector<Particle *> particles) {

    const double a = 1.0; //For N=2 the spin is always anti-parallell?
    const double omega = m_parameters[0]; const double alpha = m_parameters[1];
    const double beta  = m_parameters[2];

    double r1_squared  = 0.0;
    double r2_squared  = 0.0;
    double r12 = 0.0;

    for(int d = 0; d < m_system->getNumberOfDimensions(); d++) {
        double x1= particles[0]->getPosition()[d];
        double x2= particles[1]->getPosition()[d];
        r1_squared+=x1*x1;
        r2_squared+=x2*x2;
        r12+=(x1-x2)*(x1-x2);
    }

    r12 = sqrt(r12);

    if(m_jastrow) {
        return exp(-alpha*omega*0.5*(r1_squared+r2_squared))*exp(a*r12/(1.0+beta*r12));
    } else {
        return exp(-alpha*omega*0.5*(r1_squared+r2_squared));
    }

}

double SimpleQuantumDotWaveFunc::computeAlphaDerivative(std::vector<Particle *> particles) {

    const double C = 1.0;
    const double a = 1.0; //For N=2 the spin is always anti-parallell?
    const double omega = m_parameters[0]; const double alpha = m_parameters[1];
    const double beta  = m_parameters[2];

    double r1_squared  = 0.0;
    double r2_squared  = 0.0;
    double r12 = 0.0;

    for(int d = 0; d < m_system->getNumberOfDimensions(); d++) {
        double x1= particles[0]->getPosition()[d];
        double x2= particles[1]->getPosition()[d];
        r1_squared+=x1*x1;
        r2_squared+=x2*x2;
        r12+=(x1-x2)*(x1-x2);
    }

    r12 = sqrt(r12);

    return -omega*0.5*(r1_squared + r2_squared); //Possibly multiply with psi_t
}

double SimpleQuantumDotWaveFunc::computeBetaDerivative(std::vector<Particle *> particles) {
    const double C = 1.0;
    const double a = 1.0; //For N=2 the spin is always anti-parallell?
    const double omega = m_parameters[0]; const double alpha = m_parameters[1];
    const double beta  = m_parameters[2];

    double r1_squared  = 0.0;
    double r2_squared  = 0.0;
    double r12 = 0.0;

    for(int d = 0; d < m_system->getNumberOfDimensions(); d++) {
        double x1= particles[0]->getPosition()[d];
        double x2= particles[1]->getPosition()[d];
        r1_squared+=x1*x1;
        r2_squared+=x2*x2;
        r12+=(x1-x2)*(x1-x2);
    }

    r12 = sqrt(r12);

    return (-a*r12*r12)/((1.0+beta*r12)*(1.0+beta*r12)); //Possibly multiply with psi_t
}

std::vector<double> SimpleQuantumDotWaveFunc::computeGradient(std::vector<Particle *> particles) {

    std::vector<double> gradient(4);

    double x1 = particles[0]->getPosition()[0];
    double x2 = particles[0]->getPosition()[1];
    double y1 = particles[1]->getPosition()[0];
    double y2 = particles[1]->getPosition()[1];

    double omega = m_parameters[0];
    double alpha = m_parameters[1];
    double beta  = m_parameters[2];
    double a     = 1.0;
    double r12   = sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2));

    if(m_system->getJastrow()) {

        gradient[0] = -x1*alpha*omega + (a*(x1-x2))/(r12*(1+beta*r12)) - (a*beta*(x1-x2))/((1+beta*r12)*(1+beta*r12));
        gradient[1] = -x2*alpha*omega - (a*(x1-x2))/(r12*(1+beta*r12)) + (a*beta*(x1-x2))/((1+beta*r12)*(1+beta*r12));
        gradient[2] = -y1*alpha*omega + (a*(y1-y2))/(r12*(1+beta*r12)) - (a*beta*(y1-y2))/((1+beta*r12)*(1+beta*r12));
        gradient[3] = -y2*alpha*omega - (a*(y1-y2))/(r12*(1+beta*r12)) + (a*beta*(y1-y2))/((1+beta*r12)*(1+beta*r12));

    } else {
        gradient[0] = -x1*alpha*omega;
        gradient[1] = -x2*alpha*omega;
        gradient[2] = -y1*alpha*omega;
        gradient[3] = -y2*alpha*omega;
    }

    return gradient;
}
