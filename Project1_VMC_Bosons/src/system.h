#pragma once
#include <vector>
#include <sampler.h>

class System {
public:
    System(bool importanceSampling);
    bool metropolisStep                     ();
    bool metropolisStepImportanceSampling   ();
    bool metropolisStepAltImpSamp();
    void runMetropolisSteps                 (int numberOfMetropolisSteps);
    void setNumberOfParticles               (int numberOfParticles);
    void setNumberOfDimensions              (int numberOfDimensions);
    void setStepLength                      (double stepLength);
    void setEquilibrationFraction           (double equilibrationFraction);
    void setHamiltonian                     (class Hamiltonian* hamiltonian);
    void setWaveFunction                    (class WaveFunction* waveFunction);
    void setInitialState                    (class InitialState* initialState);
    void setAlpha                           (double alpha);
    void setBeta                            (double beta);
    void setOptimize(bool opt) {m_optimize = opt;}
    void setWriteToFile(bool wrtf) {m_writeToFile = wrtf;}
    class WaveFunction*                     getWaveFunction()   { return m_waveFunction; }
    class Hamiltonian*                      getHamiltonian()    { return m_hamiltonian; }
    class Sampler*                          getSampler()        { return m_sampler; }
    std::vector<class Particle*>            getParticles()      { return m_particles; }
    int getNumberOfParticles()              { return m_numberOfParticles; }
    int getNumberOfDimensions()             { return m_numberOfDimensions; }
    int getNumberOfMetropolisSteps()        { return m_numberOfMetropolisSteps; }
    double getEquilibrationFraction()       { return m_equilibrationFraction; }
    double getAlphaDerivativeEnergy();
    double getBetaDerivativeEnergy();
    double getEnergy();
    void setRank(int rank) {m_rank = rank;}
    void setSize(int size) {m_size = size;}
    int getRank() {return m_rank;}
    int getSize() {return m_size;}
    void printToTerminal();
    bool getOptimize() {return m_optimize;}
    bool getWriteToFile() {return m_writeToFile;}
    void setJastrow(bool jast) {m_jastrow = jast;}
    bool getJastrow() {return m_jastrow;}

private:
    int                             m_numberOfParticles = 0;
    int                             m_numberOfDimensions = 0;
    int                             m_numberOfMetropolisSteps = 0;
    int m_rank;
    int m_size;
    double                          m_equilibrationFraction = 0.0;
    double                          m_stepLength = 0.1;
    bool                            m_importanceSampling = false;
    class WaveFunction*             m_waveFunction = nullptr;
    class Hamiltonian*              m_hamiltonian = nullptr;
    class InitialState*             m_initialState = nullptr;
    class Sampler*                  m_sampler = nullptr;
    std::vector<class Particle*>    m_particles = std::vector<class Particle*>();
    bool m_optimize = false;
    bool m_writeToFile = false;
    bool m_jastrow = false;
};

