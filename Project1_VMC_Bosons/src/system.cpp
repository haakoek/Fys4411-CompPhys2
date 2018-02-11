#include "system.h"
#include <cassert>
#include "sampler.h"
#include "particle.h"
#include "WaveFunctions/wavefunction.h"
#include "Hamiltonians/hamiltonian.h"
#include "InitialStates/initialstate.h"
#include "Math/random.h"
#include <iostream>
#include <cmath>

using std::cout;
using std::endl;
using std::sqrt;

System::System(bool importanceSampling) {
    m_importanceSampling = importanceSampling;
}

//bool System::metropolisStepSlaterDet()

bool System::metropolisStep() {

    /* Perform the actual Metropolis step: Choose a particle at random and
     * change it's position by a random amount, and check if the step is
     * accepted by the Metropolis test (compare the wave function evaluated
     * at this new position with the one at the old position).
     */



    int random_particle  = Random::nextInt(m_numberOfParticles);  //Choose a random particle
    int random_dimension = Random::nextInt(m_numberOfDimensions); //Choose a random dimension
    double psi_old       = m_waveFunction->evaluate(m_particles); //Eval the wavefunction at current position

    /*
    std::vector<double> change(m_numberOfDimensions);


    for(int i = 0; i < m_numberOfDimensions; i++) {

        double r = Random::nextDouble()*2.0-1.0;
        m_particles[random_particle]->adjustPosition(r*m_stepLength,i); //Change position by r*m_stepLength
        change[i] = r;

    }
    */

    double change = Random::nextDouble()*2.0-1.0;
    m_particles[random_particle]->adjustPosition(change*m_stepLength,random_dimension);

    double psi_new = m_waveFunction->evaluate(m_particles); //Eval wavefunction at the suggested new position

    //Perform Metropolis test

    double w = (psi_new*psi_new)/(psi_old*psi_old); //Compute w = P(R')/P(R)
    double s = Random::nextDouble();

    if(w > 1) {
        return true;
    } else if(w >= s) {
        return true;
    } else {
        //Reset positions
        /*
        for(int i = 0; i < m_numberOfDimensions; i++) {
            m_particles[random_particle]->adjustPosition(-change[i]*m_stepLength,i);
        }
        */
        m_particles[random_particle]->adjustPosition(-change*m_stepLength,random_dimension);
        return false;
    }

}

bool System::metropolisStepAltImpSamp() {

    int random_particle = Random::nextInt(m_numberOfParticles);

    double psi_old = m_waveFunction->evaluate(m_particles);
    std::vector<double> F_old = m_waveFunction->computeGradient(m_particles);
    double constant_old = 2.0/psi_old;

    for(int i = 0; i < F_old.size(); i++) {
        F_old[i]*= constant_old;
    }


    std::vector<double> change(m_numberOfDimensions);

    for(int i = 0; i < m_numberOfDimensions; i++) {
        change[i] = 0.5*F_old[i]*m_stepLength + Random::nextGaussian(0.0,sqrt(m_stepLength));
        m_particles[random_particle]->adjustPosition(change[i],i);
    }

    double psi_new = m_waveFunction->evaluate(m_particles);
    std::vector<double> F_new = m_waveFunction->computeGradient(m_particles);
    double constant_new = 2.0/psi_new;

    for(int i = 0; i < F_new.size(); i++) {
        F_new[i]*= constant_new;
    }

    std::vector<double> oldPosition(m_numberOfDimensions);
    std::vector<double> newPosition(m_numberOfDimensions);

    for(int i = 0; i < m_numberOfDimensions; i++) {
        newPosition[i] = m_particles[random_particle]->getPosition()[i];
        oldPosition[i] = m_particles[random_particle]->getPosition()[i] - change[i];
    }

   double exponent = 0;


   for (int dim=0; dim < m_numberOfDimensions; dim++) {
       double term1 = - (oldPosition[dim] - newPosition[dim]  -  0.5*m_stepLength*F_new[dim]) *
                        (oldPosition[dim] - newPosition[dim]  -  0.5*m_stepLength*F_new[dim]);
       double term2 =   (-oldPosition[dim] + newPosition[dim] -  0.5*m_stepLength*F_old[dim]) *
                        (-oldPosition[dim] + newPosition[dim] -  0.5*m_stepLength*F_old[dim]);
       exponent += term1 + term2;
   }

   double greensRatio = exp(exponent / 2*m_stepLength);

   double ratio = (psi_new*psi_new)/(psi_old*psi_old)*greensRatio;

   if(ratio > 1) {
       return true;
   } else if(ratio >= Random::nextDouble()) {
       return true;
   } else {

       for(int i = 0; i < m_numberOfDimensions; i++) {
           m_particles[random_particle]->adjustPosition(-change[i],i);
       }

       return false;
   }
}

bool System::metropolisStepImportanceSampling() {

    // 1. Compute quantities in this position.

    double psi_old = m_waveFunction->evaluate(m_particles);
    std::vector<double> F_old = m_waveFunction->computeGradient(m_particles);
    double constant_old = 2.0/psi_old;
    int dim = m_numberOfDimensions;

    for(int i = 0; i < F_old.size(); i++) {
        F_old[i]*= constant_old;
    }

    // 1. end

    // 2. Propose new position

    std::vector<double> change(dim*m_numberOfParticles);

    for(int p = 0; p < m_numberOfParticles; p++) {
        for(int d = 0; d < m_numberOfDimensions; d++) {
            double dr = 0.5*m_stepLength*F_old[dim*p+d] + Random::nextGaussian(0.0,sqrt(m_stepLength));
            change[dim*p+d] = dr;
            m_particles[p]->adjustPosition(dr,d); // x + dr
        }
    }

    // 2. end

    // 3. Compute quantities in proposed position

    double psi_new = m_waveFunction->evaluate(m_particles);
    std::vector<double> F_new = m_waveFunction->computeGradient(m_particles);
    double constant_new = 2.0/psi_new;

    for(int i = 0; i < F_new.size(); i++) {
        F_new[i]*= constant_new;
    }

    // 3. end

    // 4. Metropolis Step

    //double C   = 1.0/(4.0*M_PI*0.5*m_stepLength);
    //C = pow(C,1.5*m_numberOfParticles);

    double Gyx = 0.0; // G(y,x,dt) = C*exp(-((R' - R - D*dt*F_old)^2)/(4*D*dt))
                      // v_yx = R' - R - D*dt*F_old
                      // G(y,x,dt) = const*exp(-v_yx*v_yx/(4*D*dt)), * is dot.

    double Gxy = 0.0; // G(x,y,dt) = C*exp(-((R-R'-D*dt*F_new)^2/(4*D*dt))
                      // v_xy      = R-R'-D*dt*F_new
                      // G(x,y,dt) = C*exp(-v_xy*v_xy/(4*D*dt))

    //First compute dot v_yx*v_yx = v_yx2 and v_xy*v_xy = v_xy2

    double v_yx2 = 0.0;
    double v_xy2 = 0.0;

    for(int p = 0; p < m_numberOfParticles; p++) {
        for(int d = 0; d < m_numberOfDimensions; d++) {

            double xd_n = m_particles[p]->getPosition()[d] - change[dim*p+d]; // dth-component of particle n (old)
            double yd_n = m_particles[p]->getPosition()[d];                 // dth-component of particle n (proposed)

            double Fyx_k  = F_old[dim*p+d];
            double Fxy_k  = F_new[dim*p+d];

            double v_yx_k   = yd_n - xd_n - 0.5*m_stepLength*Fyx_k;              // kth-component of v
            double v_xy_k   = xd_n - yd_n - 0.5*m_stepLength*Fxy_k;

            v_yx2 += v_yx_k*v_yx_k;
            v_xy2 += v_xy_k*v_xy_k;


        }
    }

    //cout << v_yx2 << endl;
    //cout << v_xy2 << endl;

    Gyx = exp(-v_yx2/(2*m_stepLength));
    Gxy = exp(-v_xy2/(2*m_stepLength));

    //Perform Metropolis test.

    double w = (Gxy/Gyx)*((psi_new*psi_new)/(psi_old*psi_old)); //Compute w = (G(y,x)/G(x,y))*P(R')/P(R)
    double s = Random::nextDouble();

    if(w > 1) {
        return true;
    } else if(w >= s) {
        return true;
    } else {
        //Reset positions
        for(int p = 0; p < m_numberOfParticles; p++) {
            for(int d = 0; d < m_numberOfDimensions; d++) {
                m_particles[p]->adjustPosition(-change[dim*p+d],d);
            }
        }

        return false;
    }


}

void System::runMetropolisSteps(int numberOfMetropolisSteps) {

    m_particles                 = m_initialState->getParticles();
    m_sampler                   = new Sampler(this);
    m_numberOfMetropolisSteps   = numberOfMetropolisSteps;
    m_sampler->setNumberOfMetropolisSteps(numberOfMetropolisSteps);
    bool acceptedStep = false;

    for (int i=0; i < numberOfMetropolisSteps; i++) {

        if(m_importanceSampling) {
            acceptedStep = metropolisStepAltImpSamp();
        } else {
            acceptedStep = metropolisStep();
        }

        if (i/((double) numberOfMetropolisSteps) >= m_equilibrationFraction) {
            m_sampler->sample(acceptedStep);
        }

        if(m_rank == 0) {
            if (!(i%1000)) {
                cout << "Progress: " << i/((double) numberOfMetropolisSteps) * 100 << " % \r";
                fflush(stdout);
            }
        }


    }

    m_sampler->computeAverages();
}

void System::setNumberOfParticles(int numberOfParticles) {
    m_numberOfParticles = numberOfParticles;
}

void System::setNumberOfDimensions(int numberOfDimensions) {
    m_numberOfDimensions = numberOfDimensions;
}

void System::setStepLength(double stepLength) {
    assert(stepLength >= 0);
    m_stepLength = stepLength;
}

void System::setEquilibrationFraction(double equilibrationFraction) {
    assert(equilibrationFraction >= 0);
    m_equilibrationFraction = equilibrationFraction;
}

void System::setHamiltonian(Hamiltonian* hamiltonian) {
    m_hamiltonian = hamiltonian;
}

void System::setWaveFunction(WaveFunction* waveFunction) {
    m_waveFunction = waveFunction;
}

void System::setInitialState(InitialState* initialState) {
    m_initialState = initialState;
}

void System::setAlpha(double alpha) {
    m_waveFunction->setAlpha(alpha);
}

void System::setBeta(double beta) {
    m_waveFunction->setBeta(beta);
}

double System::getAlphaDerivativeEnergy() {
    return m_sampler->getDEDalpha();
}

double System::getBetaDerivativeEnergy() {
    return m_sampler->getDEDbeta();
}

double System::getEnergy() {
    return m_sampler->getEnergy();
}

void System::printToTerminal() {
    m_sampler->printOutputToTerminal();
}


