#pragma once
#include <vector>


class WaveFunction {
public:
    WaveFunction(class System* system);
    int     getNumberOfParameters() { return m_numberOfParameters; }
    std::vector<double> getParameters() { return m_parameters; }
    virtual double evaluate(std::vector<class Particle*> particles) = 0;
    //virtual double computeDoubleDerivative(std::vector<class Particle*> particles) {};
    virtual std::vector<double> computeGradient(std::vector<Particle*> particles);
    void setAlpha(double alpha);
    void setBeta(double beta);
    void setJastrow(bool jastrow) {m_jastrow = jastrow;}
    virtual double computeAlphaDerivative(std::vector<class Particle*> particles) = 0;
    virtual double computeBetaDerivative(std::vector<class Particle*> particles);

protected:
    int     m_numberOfParameters = 0;
    std::vector<double> m_parameters = std::vector<double>();
    class System* m_system = nullptr;
    bool m_jastrow = false;
};

