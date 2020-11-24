#ifndef EVOLUTION_H
#define EVOLUTION_H
class Sim;

class Evolution{
    private:
        fftw_complex *fftwin, *fftwout;
        fftw_plan forward, backward;

        std::complex<double> *temp1;
        std::complex<double> *temp2;
        void evMomentum(std::complex<double> *in, std::complex<double> *out, double t);
	    void evSpace(std::complex<double> *in, std::complex<double> *out, double t);
  	    void evSpace(std::complex<double> *in1, std::complex<double> *in2, std::complex<double> *out, double t);
        void timeStepG(double t);
        void timeStepE(double t);

        Sim *sim;

    public:
        Evolution(Sim *_sim);
        void timeStep(double t);
        
};

#endif