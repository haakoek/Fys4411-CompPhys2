#pragma once

class Examples
{
public:
    static int nonInteractingHObruteForce(int nrOfdims, int nrOfpart, int nrOfsteps, double om);
    static int nonInteractingHOimportanceSampling(int nrOfdims, int nrOfparts, int nrOfsteps);
    static int interactingEllipticOscillator(int nrOfdims, int nrOfparts, int nrOfsteps);
    static int HeliumAtom();
    static int nonInteractingHOSteepestDescent(double initial_alpha);
    static int interactingEllipticOscillatorSteepestDescent();
    static int SimpleQuantumDot(int argc, char* argv[]);
};

