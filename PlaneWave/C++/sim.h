#ifndef SIM_H
#define SIM_H
#include <complex>
#include <fftw3.h>
#include <fstream>

//const static std::complex<double> I(std::complex<double>(0,1));

class Sim{
    private:
        int Nx, Nz;
        double xmin, xmax;
        double zmin, zmax;
        double x0, z0;
        double q0, p0;
        double a;
        std::complex<double> **Phi;
        std::complex<double> **PhiMomentum;
        std::complex<double> **PhiFuture;
        double **Phi2;

        double *x, *z;
        double *q, *p;
        //double **X, **Z;
        double dz, dx;
        double dq, dp;

        double *qshift, *pshift;

        fftw_complex *in, *out;
        fftw_plan forward, backward;

        void initMatrices();
        void initSpace();
        void initMomentum();
        void planFFT(); 
        void evMomentum(std::complex<double> *in, double t);
        void timeStep();
        
    public:
       Sim(double _xmin, double _xmax, double _zmin, 
           double _zmax, double _x0, double _z0, 
           double _q0, double _p0, double _a, 
           double _Nx, double _Nz);
       void writeWavePacket (std::ofstream& file);
       void Evolution();
};


#endif