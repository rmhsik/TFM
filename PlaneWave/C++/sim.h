#ifndef SIM_H
#define SIM_H
#include <complex>
#include <fftw3.h>
#include <fstream>
class Evolution;
//const static std::complex<double> I(std::complex<double>(0,1));

class Sim{
    private:
        std::complex<double> **PhiT;
        double **Phi2G;
        double **Phi2E;
        double **Phi2T;

        void initMatrices();
        void initSpaceTime();
        void initMomentum();
        void planFFT(); 
        
        int Nx, Nz, Nt, Nsample;
        double xmin, xmax;
        double zmin, zmax;
        double tmin, tmax;
        double x0, z0;
        double q0, p0;
        double a, m;

        double *x, *z;
        double *q, *p;
        double *t;
        double dz, dx;
        double dq, dp;
        double dt;
        double *qshift, *pshift;

        fftw_complex *in, *out;
        fftw_plan forward, backward;

        std::complex<double> **PhiG;
        std::complex<double> **PhiMomentumG;
        std::complex<double> **PhiFutureG;
        std::complex<double> **PhiE;
        std::complex<double> **PhiMomentumE;
        std::complex<double> **PhiFutureE;

        Evolution *evOperator;

    protected:
        
    public:
       Sim(double _xmin, double _xmax, double _zmin, 
           double _zmax, double _tmin, double _tmax,
           double _x0, double _z0, 
           double _q0, double _p0, double _a, 
           int _Nx, int _Nz, int _Nt, int _Nsample);
       void writeWavePacket (int j);
       void Benchmark();
       void Run();
       void write2File(std::complex<double> **Phi, double **Phi2,
                       const char *name, int j);
        void write2FileT(int j);
       friend class Evolution;
};


#endif
