#include "examples.h"

int main(int argc, char* argv[]) {
    //return Examples::nonInteractingHObruteForce(1,10,1e5,1.0);
    //return Examples::nonInteractingHOimportanceSampling(1,10,1e5);
    //return Examples::HeliumAtom();
    return Examples::interactingEllipticOscillator(3,10,1e4);
    //return Examples::nonInteractingHOSteepestDescent(0.7);
    //return Examples::interactingEllipticOscillatorSteepestDescent();
    //return Examples::SimpleQuantumDot(argc,argv);
}
