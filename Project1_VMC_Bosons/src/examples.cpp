#include "examples.h"
#include <iostream>
#include <iomanip>
#include "system.h"
#include "particle.h"
#include "WaveFunctions/wavefunction.h"
#include "WaveFunctions/simplegaussian.h"
#include "WaveFunctions/interactinggaussian.h"
#include "WaveFunctions/heliumwavefunction.h"
#include "WaveFunctions/simplequantumdotwavefunc.h"
#include "Hamiltonians/heliumhamiltonian.h"
#include "Hamiltonians/hamiltonian.h"
#include "Hamiltonians/harmonicoscillator.h"
#include "InitialStates/initialstate.h"
#include "InitialStates/randomuniform.h"
#include "Math/random.h"
#include <ctime>
#include "Hamiltonians/ellipticinteractingoscillator.h"
#include "Hamiltonians/simplequantumdothamiltonian.h"
#include "steepestdescent.h"
#include <mpi/mpi.h>
#include <cstdlib>

using namespace std;

int Examples::nonInteractingHObruteForce(int nrOfdims, int nrOfpart,int nrOfsteps, double om) {

    int numberOfDimensions  = nrOfdims;
    int numberOfParticles   = nrOfpart;
    int numberOfSteps       = (int) nrOfsteps;
    double omega            = om;                    // Oscillator frequency.
    double alpha            = 0.5*omega;            // Variational parameter.
    double stepLength       = 2.5;                 // Metropolis step length.
    double equilibration    = 0.1;                // Amount of the total steps used for equilibration.
    bool analytical         = true;             // Decide whether to use numerical diff or analytical expression for local energy
    bool imp_sampling       = false;

    System* system =                     new System(imp_sampling);
    system->setHamiltonian              (new HarmonicOscillator(system, omega,analytical));
    system->setWaveFunction             (new SimpleGaussian(system, alpha));
    system->setInitialState             (new RandomUniform(system, numberOfDimensions, numberOfParticles));


    system->setEquilibrationFraction    (equilibration);
    system->setStepLength               (stepLength);

    clock_t begin = clock();

    system->runMetropolisSteps(numberOfSteps);

    clock_t end   = clock();

    cout << "Execution time: " << double (end - begin) / CLOCKS_PER_SEC << " seconds";

    system->printToTerminal();

    return 0;
}

int Examples::nonInteractingHOimportanceSampling(int nrOfdims, int nrOfparts, int nrOfsteps) {

    int numberOfDimensions  = nrOfdims;
    int numberOfParticles   = nrOfparts;
    int numberOfSteps       = (int) nrOfsteps;
    double omega            = 1.0;          // Oscillator frequency.
    double alpha            = 0.5;          // Variational parameter.
    double stepLength       = 0.1;          // Metropolis step length.
    double equilibration    = 0.1;          // Amount of the total steps used for equilibration.
    bool analytical         = false;         // Decide whether to use numerical diff or analytical expression for local energy
    bool imp_sampling       = true;

    System* system =                     new System(imp_sampling);
    system->setHamiltonian              (new HarmonicOscillator(system, omega,analytical));
    system->setWaveFunction             (new SimpleGaussian(system, alpha));
    system->setInitialState             (new RandomUniform(system, numberOfDimensions, numberOfParticles));


    system->setEquilibrationFraction    (equilibration);
    system->setStepLength               (stepLength);

    clock_t begin = clock();
    system->runMetropolisSteps          (numberOfSteps);
    clock_t end = clock();

    double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;

    cout << "Execution time: " << elapsed_secs << " seconds." << endl;
    if(analytical) {
        cout << "Analytical local energy " << endl;
    } else {
        cout << "Numerical local energy " << endl;
    }

    system->printToTerminal();

    return 0;

}

int Examples::HeliumAtom() {
    int numberOfDimensions = 3;
    int numberOfParticles  = 2;
    int numberOfSteps      = 1e6;
    double alpha           = 1.843;
    double beta            = 0.347;
    double stepLength      = 1.5;
    double equilibration   = 0.1;
    bool imp_sampling      = false;

    System* system =                     new System(imp_sampling);
    system->setHamiltonian              (new HeliumHamiltonian(system,true));
    system->setWaveFunction             (new HeliumWaveFunction(system,alpha,beta));

    system->setInitialState             (new RandomUniform(system, numberOfDimensions, numberOfParticles));
    system->setEquilibrationFraction    (equilibration);
    system->setStepLength               (stepLength);

    clock_t begin = clock();

    system->runMetropolisSteps(numberOfSteps);

    clock_t end   = clock();

    cout << "Execution time: " << double (end - begin) / CLOCKS_PER_SEC << " seconds";

    system->printToTerminal();

    return 0;
}


int Examples::interactingEllipticOscillator(int nrOfdims, int nrOfparts, int nrOfsteps) {

    int numberOfDimensions  = nrOfdims;
    int numberOfParticles   = nrOfparts;
    int numberOfSteps       = (int) nrOfsteps;
    double gamma            = 2.82843;
    double beta             = gamma;
    double alpha            = 0.5;          // Variational parameter.
    double stepLength       = 1.5;          // Metropolis step length.
    double equilibration    = 0.2;          // Amount of the total steps used for equilibration.
    bool imp_sampling       = false;

    System* system =                     new System(imp_sampling);
    system->setHamiltonian              (new EllipticInteractingOscillator(system,gamma));
    system->setWaveFunction             (new InteractingGaussian(system,alpha,beta));
    system->setInitialState             (new RandomUniform(system, numberOfDimensions, numberOfParticles));

    system->setEquilibrationFraction    (equilibration);
    system->setStepLength               (stepLength);

    clock_t begin = clock();

    system->runMetropolisSteps(numberOfSteps);

    clock_t end   = clock();

    cout << "Execution time: " << double (end - begin) / CLOCKS_PER_SEC << " seconds";

    system->printToTerminal();

    return 0;

}

int Examples::nonInteractingHOSteepestDescent(double initial_alpha) {
    int numberOfDimensions  = 3;
    int numberOfParticles   = 10;
    double omega            = 1.0;                     // Oscillator frequency.
    double alpha0           = initial_alpha;          // Variational parameter.
    double stepLength       = 2.0;                   // Metropolis step length.
    double equilibration    = 0.1;                  // Amount of the total steps used for equilibration.
    bool analytical         = true;                // Decide whether to use numerical diff or analytical expression for local energy
    bool imp_sampling       = false;

    System* system =                     new System(imp_sampling);
    system->setHamiltonian              (new HarmonicOscillator(system, omega,analytical));
    system->setWaveFunction             (new SimpleGaussian(system, alpha0));
    system->setInitialState             (new RandomUniform(system, numberOfDimensions, numberOfParticles));

    system->setEquilibrationFraction    (equilibration);
    system->setStepLength               (stepLength);

    SteepestDescent* sd = new SteepestDescent(system,1e4,numberOfParticles,numberOfDimensions);
    double optimal_alpha = sd->findOptimalAlpha(alpha0);


    int nr_of_iterations = sd->getNrOfIters();
    cout << "Number iterations before optimal alpha found: " << nr_of_iterations << endl;
    cout << "Optimal alpha: " << optimal_alpha << endl;

    return 0;
}

int Examples::interactingEllipticOscillatorSteepestDescent() {
    int numberOfDimensions  = 3;
    int numberOfParticles   = 50;
    int numberOfSteps       = (int) 1e3;
    double gamma            = 2.82843;
    double beta             = gamma;
    double alpha            = 0.52;          // Variational parameter.
    double stepLength       = 1.5;          // Metropolis step length.
    double equilibration    = 0.2;          // Amount of the total steps used for equilibration.
    bool imp_sampling       = false;

    System* system =                     new System(imp_sampling);
    system->setHamiltonian              (new EllipticInteractingOscillator(system,gamma));
    system->setWaveFunction             (new InteractingGaussian(system,alpha,beta));
    system->setInitialState             (new RandomUniform(system, numberOfDimensions, numberOfParticles));

    system->setEquilibrationFraction    (equilibration);
    system->setStepLength               (stepLength);

    SteepestDescent* sd = new SteepestDescent(system,numberOfSteps,numberOfParticles,numberOfDimensions);
    double optimal_alpha = sd->findOptimalAlpha(alpha);

    cout << optimal_alpha << endl;

    return 0;
}

int Examples::SimpleQuantumDot(int argc, char* argv[]) {

    MPI_Init (&argc, &argv);	/* starts MPI */
    int rank, size;
    MPI_Comm_rank (MPI_COMM_WORLD, &rank);	/* get current process id */
    MPI_Comm_size (MPI_COMM_WORLD, &size);	/* get number of processes */

    Random::setSeed((long) -rank*10);
    cout << "rank;: " << rank << ", seed:" << Random::getSeed() << " first=" << Random::nextDouble() << endl;
    int numberOfDimensions = 2;
    int numberOfParticles  = 2;
    int numberOfSteps = (int) 1e6;
    double stepLength = 0.0001;
    double equilibration = 0.1;
    bool imp_sampling = true;
    bool optimize     = false;
    bool writeToFile  = true;

    double omega, alpha, beta;
    int int_jastrow;
    bool jastrow;

    if(argc > 1) {
        omega = atof(argv[1]);
        alpha = atof(argv[2]);
        beta  = atof(argv[3]);
        int_jastrow = atoi(argv[4]);

        if(int_jastrow == 1) {
            jastrow = true;
        } else {
            jastrow = false;
        }

    }

    System* system = new System(imp_sampling);
    system->setRank(rank);
    system->setSize(size);
    system->setHamiltonian(new SimpleQuantumDotHamiltonian(system,true));
    system->setWaveFunction(new SimpleQuantumDotWaveFunc(system,omega,alpha,beta));
    system->setInitialState(new RandomUniform(system,numberOfDimensions,numberOfParticles));
    system->setEquilibrationFraction(equilibration);
    system->setStepLength(stepLength);
    system->setOptimize(optimize);
    system->setWriteToFile(writeToFile);
    system->setJastrow(jastrow);
    system->getWaveFunction()->setJastrow(jastrow);

    if(optimize) {

        SteepestDescent* sd = new SteepestDescent(system,1e6,numberOfParticles,numberOfDimensions);

        if(jastrow) {
            std::vector<double> parameters(2);
            parameters[0] = alpha;
            parameters[1] = beta;
            sd->altOptimize(parameters);
        } else {
            sd->findOptimalAlpha(alpha);
        }

    } else {


        clock_t begin = clock();

        system->runMetropolisSteps(numberOfSteps);

        clock_t end   = clock();

        system->printToTerminal();

        cout << "Execution time: " << double (end - begin) / CLOCKS_PER_SEC << " seconds" << endl;
    }

    MPI_Finalize();

    return 0;

}

