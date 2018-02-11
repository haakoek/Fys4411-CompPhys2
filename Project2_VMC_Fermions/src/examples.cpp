#include "examples.h"
#include <iostream>
#include <iomanip>
#include "system.h"
#include "particle.h"
#include "WaveFunctions/wavefunction.h"
#include "Hamiltonians/hamiltonian.h"
#include "InitialStates/initialstate.h"
#include "InitialStates/randomuniform.h"
#include "Math/random.h"
#include <ctime>
#include "steepestdescent.h"
#include <WaveFunctions/slaterwavefunction.h>
#include <WaveFunctions/wavefunction.h>
#include <SingleParticleWaveFunctions/singleparticleharmonicoscillator.h>
#include <SingleParticleWaveFunctions/singleparticlewavefunctions.h>
#include "Hamiltonians/simplequantumdothamiltonian.h"

using namespace std;

int Examples::TwoParticles() {

    int numberOfDimensions = 2;
    int numberOfParticles  = 2;
    int numberOfSteps = (int) 1e6;
    double omega = 0.01;
    double alpha            = 0.946275;          // variational parameter 1
    double beta             = 0.300194;        // variational parameter 2
    double stepLength = 1.7;
    double timeStep   = 0.001;
    double equilibration = 0.1;

    bool interaction  = true;
    bool with_jastrow = true;
    bool imp_sampling = true;
    bool writeToFile  = false;

    System* system =                     new System(imp_sampling, with_jastrow,writeToFile);
    system->setHamiltonian(new SimpleQuantumDotHamiltonian(system,with_jastrow,interaction));

    system->setInitialState             (new RandomUniform(system, numberOfDimensions, numberOfParticles));
    WaveFunction* wf = new SlaterWaveFunction(system,omega,alpha,beta,new SingleParticleHarmonicOscillator(omega,alpha));
    system->setWaveFunction(wf);

    system->setEquilibrationFraction(equilibration);
    system->setStepLength(stepLength);
    system->setTimeStepImportanceSampling(timeStep);


    std::vector<double> parameters(2);
    parameters[0] = alpha;
    parameters[1] = beta;
    SteepestDescent* sd = new SteepestDescent(system,1e6,numberOfParticles,numberOfDimensions);

    //sd->optimize(parameters);
    sd->altOptimize(parameters);

    system->setAlpha(sd->getOptAlpha());
    system->setBeta(sd->getOptBeta());

    //system->setWriteToFile(true);

    /*
    clock_t begin = clock();

    system->runMetropolisSteps(numberOfSteps);

    clock_t end   = clock();

    system->printToTerminal();
    cout << "Execution time: " << double (end - begin) / CLOCKS_PER_SEC << " seconds" << endl;
    */

    return 0;
}

int Examples::SixParticles() {

    int numberOfDimensions = 2;
    int numberOfParticles  = 6;
    int numberOfSteps = (int) 1e6;
    double omega = 0.01;
    double alpha            = 0.849742;          // variational parameter 1
    double beta             = 0.155566;        // variational parameter 2
    double stepLength = 1.7;
    double timeStep   = 0.001;
    double equilibration = 0.1;

    bool interaction  = true;
    bool with_jastrow = true;
    bool imp_sampling = true;
    bool writeToFile  = false;

    System* system =                     new System(imp_sampling, with_jastrow,writeToFile);
    system->setHamiltonian(new SimpleQuantumDotHamiltonian(system,with_jastrow,interaction));

    system->setInitialState             (new RandomUniform(system, numberOfDimensions, numberOfParticles));
    WaveFunction* wf = new SlaterWaveFunction(system,omega,alpha,beta,new SingleParticleHarmonicOscillator(omega,alpha));
    system->setWaveFunction(wf);

    system->setEquilibrationFraction(equilibration);
    system->setStepLength(stepLength);
    system->setTimeStepImportanceSampling(timeStep);


    std::vector<double> parameters(2);
    parameters[0] = alpha;
    parameters[1] = beta;
    SteepestDescent* sd = new SteepestDescent(system,1e6,numberOfParticles,numberOfDimensions);

    //sd->optimize(parameters);
    sd->altOptimize(parameters);

    system->setAlpha(sd->getOptAlpha());
    system->setBeta(sd->getOptBeta());

    /*
    clock_t begin = clock();

    system->runMetropolisSteps(numberOfSteps);

    clock_t end   = clock();

    system->printToTerminal();
    cout << "Execution time: " << double (end - begin) / CLOCKS_PER_SEC << " seconds" << endl;
    */

    return 0;

}

int Examples::TwelveParticles() {


    int numberOfDimensions = 2;
    int numberOfParticles  = 12;
    int numberOfSteps = (int) 1e6;
    double omega = 0.01;
    double alpha            = 0.84;          // variational parameter 1
    double beta             = 0.10;        // variational parameter 2
    double stepLength = 1.7;
    double timeStep   = 0.001;
    double equilibration = 0.1;

    bool interaction  = true;
    bool with_jastrow = true;
    bool imp_sampling = true;
    bool writeToFile  = false;

    System* system =                     new System(imp_sampling, with_jastrow,writeToFile);
    system->setHamiltonian(new SimpleQuantumDotHamiltonian(system,with_jastrow,interaction));

    system->setInitialState             (new RandomUniform(system, numberOfDimensions, numberOfParticles));
    WaveFunction* wf = new SlaterWaveFunction(system,omega,alpha,beta,new SingleParticleHarmonicOscillator(omega,alpha));
    system->setWaveFunction(wf);

    system->setEquilibrationFraction(equilibration);
    system->setStepLength(stepLength);
    system->setTimeStepImportanceSampling(timeStep);


    std::vector<double> parameters(2);
    parameters[0] = alpha;
    parameters[1] = beta;
    SteepestDescent* sd = new SteepestDescent(system,1e6,numberOfParticles,numberOfDimensions);

    //sd->optimize(parameters);
    sd->altOptimize(parameters);

    system->setAlpha(sd->getOptAlpha());
    system->setBeta(sd->getOptBeta());

    /*
    clock_t begin = clock();

    system->runMetropolisSteps(numberOfSteps);

    clock_t end   = clock();

    system->printToTerminal();
    cout << "Execution time: " << double (end - begin) / CLOCKS_PER_SEC << " seconds" << endl;
    */

    return 0;


}

int Examples::TwentyParticles() {

    int numberOfDimensions = 2;
    int numberOfParticles  = 20;
    int numberOfSteps = (int) 1e5;
    double omega = 0.05;
    double alpha            = 0.842;          // variational parameter 1
    double beta             = 0.153;        // variational parameter 2
    double stepLength = 1.7;
    double timeStep   = 0.001;
    double equilibration = 0.1;

    bool interaction  = true;
    bool with_jastrow = true;
    bool imp_sampling = true;
    bool writeToFile  = false;

    System* system =                     new System(imp_sampling, with_jastrow,writeToFile);
    system->setHamiltonian(new SimpleQuantumDotHamiltonian(system,with_jastrow,interaction));

    system->setInitialState             (new RandomUniform(system, numberOfDimensions, numberOfParticles));
    WaveFunction* wf = new SlaterWaveFunction(system,omega,alpha,beta,new SingleParticleHarmonicOscillator(omega,alpha));
    system->setWaveFunction(wf);

    system->setEquilibrationFraction(equilibration);
    system->setStepLength(stepLength);
    system->setTimeStepImportanceSampling(timeStep);


    std::vector<double> parameters(2);
    parameters[0] = alpha;
    parameters[1] = beta;
    SteepestDescent* sd = new SteepestDescent(system,1e6,numberOfParticles,numberOfDimensions);

    //sd->optimize(parameters);
    sd->altOptimize(parameters);

    system->setAlpha(sd->getOptAlpha());
    system->setBeta(sd->getOptBeta());

    /*
    clock_t begin = clock();

    system->runMetropolisSteps(numberOfSteps);

    clock_t end   = clock();

    system->printToTerminal();
    cout << "Execution time: " << double (end - begin) / CLOCKS_PER_SEC << " seconds" << endl;
    */

    return 0;

}
