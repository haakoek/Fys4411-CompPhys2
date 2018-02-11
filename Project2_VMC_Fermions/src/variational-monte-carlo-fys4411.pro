TEMPLATE = app
CONFIG  += console c++11
CONFIG  -= app_bundle
CONFIG  -= qt
QMAKE_CXXFLAGS_RELEASE -= -O2
QMAKE_CXXFLAGS_RELEASE += -O3


SOURCES += main.cpp \
    system.cpp \
    Hamiltonians/hamiltonian.cpp \
    particle.cpp \
    WaveFunctions/wavefunction.cpp \
    InitialStates/initialstate.cpp \
    InitialStates/randomuniform.cpp \
    Math/random.cpp \
    sampler.cpp \
    examples.cpp \
    steepestdescent.cpp \
    SingleParticleWaveFunctions/singleparticlewavefunctions.cpp \
    SingleParticleWaveFunctions/singleparticleharmonicoscillator.cpp \
    Math/hermitepolynomials.cpp \
    WaveFunctions/slaterwavefunction.cpp \
    Hamiltonians/simplequantumdothamiltonian.cpp

HEADERS += \
    system.h \
    Hamiltonians/hamiltonian.h \
    particle.h \
    WaveFunctions/wavefunction.h \
    InitialStates/initialstate.h \
    InitialStates/randomuniform.h \
    Math/random.h \
    sampler.h \
    examples.h \
    steepestdescent.h \
    SingleParticleWaveFunctions/singleparticlewavefunctions.h \
    SingleParticleWaveFunctions/singleparticleharmonicoscillator.h \
    Math/hermitepolynomials.h \
    WaveFunctions/slaterwavefunction.h \
    Hamiltonians/simplequantumdothamiltonian.h

LIBS += -llapack -lblas -larmadillo

