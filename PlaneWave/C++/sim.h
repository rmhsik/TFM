#ifndef SIM_H
#define SIM_H
#include <complex>
#include <fftw3.h>
#include <fstream>

//const static std::complex<double> I(std::complex<double>(0,1));

class Sim{
    private:
        int Nx, Nz, Nt, Nsample;
        double xmin, xmax;
        double zmin, zmax;
        double tmin, tmax;
        double x0, z0;
        double q0, p0;
        double a;
        std::complex<double> **Phi;
        std::complex<double> **PhiMomentum;
        std::complex<double> **PhiFuture;
        double **Phi2;
        double ***Data;
        
        double *x, *z;
        double *q, *p;
        double *t;
        double dz, dx;
        double dq, dp;
        double dt;
        double *qshift, *pshift;

        fftw_complex *in, *out;
        fftw_plan forward, backward;

        void initMatrices();
        void initSpaceTime();
        void initMomentum();
        void planFFT(); 
        void evMomentum(std::complex<double> *in, double t);
        void timeStep(double t);
        
    public:
       Sim(double _xmin, double _xmax, double _zmin, 
           double _zmax, double _tmin, double _tmax,
           double _x0, double _z0, 
           double _q0, double _p0, double _a, 
           int _Nx, int _Nz, int _Nt, int _Nsample);
       void writeWavePacket (int j);
       void Benchmark();
       void Evolution();
};


#endif