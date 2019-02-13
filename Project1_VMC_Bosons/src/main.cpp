#include "examples.h"

int main(int argc, char* argv[]) {

    //alpha = 0.5*omega, thus omega = 1 -> alpha = 0.5
    double alpha = 0.5;
    double omega = 1.0;
    return Examples::nonInteractingHObruteForce(3,10,1e5,omega,alpha);

    //return Examples::nonInteractingHOimportanceSampling(1,10,1e5);
    //return Examples::HeliumAtom();
    //return Examples::interactingEllipticOscillator(3,10,1e4);
    //return Examples::nonInteractingHOSteepestDescent(0.7);
    //return Examples::interactingEllipticOscillatorSteepestDescent();
    //return Examples::SimpleQuantumDot(argc,argv);
}
