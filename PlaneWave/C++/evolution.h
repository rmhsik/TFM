#ifndef EVOLUTION_H
#define EVOLUTION_H
#include "sim.h"

class Evolution{
    private:
        fftw_complex *fftwin, *fftwout;
        fftw_plan forward, backward;
        void evMomentum(std::complex<double> *in, double t);
	    void evSpace(std::complex<double> *in, double t);
        void timeStepG(double t);
        void timeStepE(double t);
        Sim *sim;

    public:
        Evolution(Sim *_sim);
        void timeStep(double t);
        
};

#endif