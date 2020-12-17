#ifndef EVOLUTION_H
#define EVOLUTION_H
class Sim;

class Evolution{
    private:
        fftw_complex *fftwin, *fftwout;
        fftw_plan forward, backward;
        double omega, l;
        std::complex<double> *temp1;
        std::complex<double> *temp2;
        double envelope(double x, double z);
        double potential(double x, double z, double t);
        void evMomentumG(std::complex<double> *in, std::complex<double> *out, double t);
        void evMomentumE(std::complex<double> *in, std::complex<double> *out, double t);
  	    void evSpaceG(std::complex<double> *in1, std::complex<double> *in2, std::complex<double> *out, double t);
  	    void evSpaceE(std::complex<double> *in1, std::complex<double> *in2, std::complex<double> *out, double t);
        void timeStepG(double t);
        void timeStepE(double t);

        Sim *sim;

    public:
        Evolution(Sim *_sim);
        void timeStep(double t);
        
};

#endif