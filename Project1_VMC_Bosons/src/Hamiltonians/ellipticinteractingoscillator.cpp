#include "ellipticinteractingoscillator.h"
#include <cassert>
#include <iostream>
#include "../system.h"
#include "../particle.h"
#include "../WaveFunctions/wavefunction.h"
#include <cmath>

using namespace std;

EllipticInteractingOscillator::EllipticInteractingOscillator(System* system, double gamma) :
    Hamiltonian(system) {
    assert(gamma > 0);
    m_gamma = gamma;
}

double EllipticInteractingOscillator::computeLocalEnergy(std::vector<Particle *> particles) {
    double kineticEnergy   = computeKineticEnergy(particles);
    double potentialEnergy = computePotentialEnergy(particles);
    return kineticEnergy + potentialEnergy;
}

double EllipticInteractingOscillator::computeKineticEnergy(std::vector<Particle *> particles) {

    const double alpha = m_system->getWaveFunction()->getParameters()[0];
    const double beta  = m_system->getWaveFunction()->getParameters()[1];
    const double a = 0.0043;

    std::vector<double> gradient = m_system->getWaveFunction()->computeGradient(particles);


    double laplacian = 0.0;

    for(int k = 0; k < m_system->getNumberOfParticles(); k++) {

        const double xk = particles[k]->getPosition()[0];
        const double yk = particles[k]->getPosition()[1];
        const double zk = particles[k]->getPosition()[2];

        const double phi_rk = exp(-alpha*(xk*xk + yk*yk + beta*zk*zk));
        const double ddphi_k = (-4.0*alpha - 2.0*alpha*beta + (2.0*alpha*xk)*(2.0*alpha*xk)
                                + (2.0*alpha*yk)*(2.0*alpha*yk)
                                + (2.0*alpha*beta*zk)*(2.0*alpha*beta*zk))*phi_rk;


        std::vector<double> sum1(3);
        double sum2 = 0.0;
        double sum3 = 0.0;

        for(int j = 0; j < m_system->getNumberOfParticles(); j++) {



            if(j != k) {

                double r_kj = 0.0;
                const double xj = particles[j]->getPosition()[0];
                const double yj = particles[j]->getPosition()[1];
                const double zj = particles[j]->getPosition()[2];

                r_kj = (xk-xj)*(xk-xj) + (yk-yj)*(yk-yj) + (zk-zj)*(zk-zj);
                r_kj = sqrt(r_kj);

                double du_kj = a/(r_kj*r_kj - a*r_kj);
                double ddu_kj = -(a*(2*r_kj - a))/((r_kj*r_kj - a*r_kj)*(r_kj*r_kj - a*r_kj));

                for(int i = 0; i < m_system->getNumberOfParticles(); i++) {
                    if(i != k) {

                        double r_ki     = 0.0;
                        const double xi = particles[i]->getPosition()[0];
                        const double yi = particles[i]->getPosition()[1];
                        const double zi = particles[i]->getPosition()[2];

                        r_ki = (xk-xi)*(xk-xi) + (yk-yi)*(yk-yi) + (zk-zi)*(zk-zi);
                        r_ki = sqrt(r_ki);

                        double du_ki = a/(r_ki*r_ki - a*r_ki);

                        double rkri_dot_rkrj = (xk-xi)*(xk-xj) + (yk-yi)*(yk-yj) + (zk-zi)*(zk-zj);

                        sum2+= ((rkri_dot_rkrj)/(r_ki*r_kj))*du_ki*du_kj;

                    }
                }


                sum1[0]+= (xk-xj)*(du_kj/r_kj);
                sum1[1]+= (yk-yj)*(du_kj/r_kj);
                sum1[2]+= (zk-zj)*(du_kj/r_kj);

                sum3   += ddu_kj + (2.0/r_kj)*du_kj;

            }
        }

        const double grad_xk = gradient[m_system->getNumberOfDimensions()*k+0];
        const double grad_yk = gradient[m_system->getNumberOfDimensions()*k+1];
        const double grad_zk = gradient[m_system->getNumberOfDimensions()*k+2];

        double grad_rk_dot_sum1 = 2.0*(grad_xk*sum1[0] + grad_yk*sum1[1] + grad_zk*sum1[2]);
        laplacian+= (ddphi_k/phi_rk) + (grad_rk_dot_sum1/phi_rk) + sum2 + sum3;

    }

    return -0.5*laplacian;
}

double EllipticInteractingOscillator::computePotentialEnergy(std::vector<Particle *> particles) {

    double potentialEnergy = 0.0;

    for(int i = 0; i < m_system->getNumberOfParticles(); i++) {
        double r2 = 0;

        for(int j=0; j < m_system->getNumberOfDimensions(); j++) {
            if(j < 2) {
                r2 += particles[i]->getPosition()[j]*particles[i]->getPosition()[j];
            } else {
                r2 += m_gamma*m_gamma*particles[i]->getPosition()[j]*particles[i]->getPosition()[j];
            }
        }

        potentialEnergy+=r2;
    }

    potentialEnergy = 0.5*potentialEnergy;


    for(int i = 0; i < m_system->getNumberOfParticles(); i++) {
        for(int j = i+1; j < m_system->getNumberOfParticles(); j++) {
            potentialEnergy+=interaction(m_system->getParticles()[i],m_system->getParticles()[j]);
        }
    }


    return potentialEnergy;
}

double EllipticInteractingOscillator::interaction(Particle* pi, Particle* pj) {

    double a = 0.0043;

    double rij = 0.0;

    for(int d = 0; d < m_system->getNumberOfDimensions(); d++) {
        rij+= (pi->getPosition()[d]-pj->getPosition()[d])*(pi->getPosition()[d]-pj->getPosition()[d]);
    }

    rij = sqrt(rij);

    if(rij <= a) {
        return 1e10;
    } else {
        return 0.0;
    }

}



