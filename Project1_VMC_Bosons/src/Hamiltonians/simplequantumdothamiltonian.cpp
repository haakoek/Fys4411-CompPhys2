#include "simplequantumdothamiltonian.h"
#include <cassert>
#include <iostream>
#include "../system.h"
#include "../particle.h"
#include "../WaveFunctions/wavefunction.h"
#include <cmath>

using namespace std;

SimpleQuantumDotHamiltonian::SimpleQuantumDotHamiltonian(System* system, bool interaction) :
    Hamiltonian(system) {
        m_interaction = interaction;

}

double SimpleQuantumDotHamiltonian::computeLocalEnergy(std::vector<Particle *> particles) {

    double r12 = 0.0;
    double r1_squared  = 0;
    double r2_squared  = 0;


    const double omega = m_system->getWaveFunction()->getParameters()[0];
    const double alpha = m_system->getWaveFunction()->getParameters()[1];
    const double beta  = m_system->getWaveFunction()->getParameters()[2];

    const double a = 1.0;
    const double A = -alpha*omega*0.5;

    const double x1 = particles[0]->getPosition()[0];
    const double y1 = particles[0]->getPosition()[1];
    const double x2 = particles[1]->getPosition()[0];
    const double y2 = particles[1]->getPosition()[1];


    r1_squared = x1*x1 + y1*y1;
    r2_squared = x2*x2 + y2*y2;
    r12        = sqrt((x2-x1)*(x2-x1) + (y2-y1)*(y2-y1));

    double dx1, dy1, dx2, dy2, ddx1, ddx2, ddy1, ddy2, term1, term2;

    if(m_system->getJastrow()) {

        dx1 = 2.0*A*x1 + a*(x1-x2)/(r12*(1.0+beta*r12)*(1.0+beta*r12));
        dy1 = 2.0*A*y1 + a*(y1-y2)/(r12*(1.0+beta*r12)*(1.0+beta*r12));
        dx2 = 2.0*A*x2 - a*(x1-x2)/(r12*(1.0+beta*r12)*(1.0+beta*r12));
        dy2 = 2.0*A*y2 - a*(y1-y2)/(r12*(1.0+beta*r12)*(1.0+beta*r12));

        term1 = a/(r12*(1.0+beta*r12)*(1.0+beta*r12)) - ((a*(x1-x2)*(x1-x2))/(r12*r12*r12*(1.0+beta*r12)*(1.0+beta*r12)*(1.0+beta*r12)))*(1.0+3.0*beta*r12);
        term2 = a/(r12*(1.0+beta*r12)*(1.0+beta*r12)) - ((a*(y1-y2)*(y1-y2))/(r12*r12*r12*(1.0+beta*r12)*(1.0+beta*r12)*(1.0+beta*r12)))*(1.0+3.0*beta*r12);

        ddx1 = 2.0*A + 2.0*A*x1*dx1 + term1 + (a*(x1-x2)/(r12*(1.0+beta*r12)*(1.0+beta*r12)))*dx1;
        ddx2 = 2.0*A + 2.0*A*x2*dx2 + term1 - (a*(x1-x2)/(r12*(1.0+beta*r12)*(1.0+beta*r12)))*dx2;
        ddy1 = 2.0*A + 2.0*A*y1*dy1 + term2 + (a*(y1-y2)/(r12*(1.0+beta*r12)*(1.0+beta*r12)))*dy1;
        ddy2 = 2.0*A + 2.0*A*y2*dy2 + term2 - (a*(y1-y2)/(r12*(1.0+beta*r12)*(1.0+beta*r12)))*dy2;

    } else {

        dx1 = 2.0*A*x1;
        dy1 = 2.0*A*y1;
        dx2 = 2.0*A*x2;
        dy2 = 2.0*A*y2;

        ddx1 = 2.0*A + 2.0*A*x1*dx1;
        ddx2 = 2.0*A + 2.0*A*x2*dx2;
        ddy1 = 2.0*A + 2.0*A*y1*dy1;
        ddy2 = 2.0*A + 2.0*A*y2*dy2;

    }

    m_potentialEnergy =  0.5*omega*omega*(r1_squared+r2_squared) + 1.0/r12;
    m_kineticEnergy   = -0.5*(ddx1+ddx2+ddy1+ddy2);

    return m_kineticEnergy + m_potentialEnergy;

}

double SimpleQuantumDotHamiltonian::getKineticEnergy(std::vector<Particle *> particles) {
    return m_kineticEnergy; //*m_system->getWaveFunction()->evaluate(particles); //Multiply with psi_t
}

double SimpleQuantumDotHamiltonian::getPotentialEnergy(std::vector<Particle *> particles) {
    return m_potentialEnergy; //*m_system->getWaveFunction()->evaluate(particles); //Multiply with psi_t
}
