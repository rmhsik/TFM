#include <iostream>
#include <cstring>
#include "sim.h"
#include "utils.h"

Sim::Sim(double _xmin, double _xmax, double _zmin,
         double _zmax, double _x0, double _z0,
          double _q0, double _p0, 
         double _a, double _Nx, double _Nz){
    
    Nx = _Nx;
    Nz = _Nz;
    xmin = _xmin;
    xmax = _xmax;
    zmin = _zmin;
    zmax = _zmax;
    x0 = _x0;
    z0 = _z0;
    q0 = _q0;
    p0 = _p0;
    a = _a;

    initSpace();
    initMomentum();
    initMatrices();
    planFFT();
}

void Sim::writeWavePacket(std::ofstream& file){
    waveSqr(Phi, Phi2, Nx, Nz);
    for (int i=0;i<Nx; i++){
        for (int j=0; j<Nz; j++){
            file << Phi2[i][j]<<" ";
        }
        file<<"\n";
    }
}

void Sim::Evolution(){
    timeStep();
}

///////////////////////////////////////////////////////////////////

void Sim::initSpace(){
    x = (double*)malloc(Nx*sizeof(double));
    z = (double*)malloc(Nz*sizeof(double));
    linspace(xmin,xmax,Nx,x,&dx);
    linspace(zmin,zmax,Nz,z,&dz);
}

void Sim::initMomentum(){
    q = (double*)malloc(Nx*sizeof(double));
    p = (double*)malloc(Nz*sizeof(double));

    linspace(-Nx/2,Nx/2,Nx,p,&dp);
    linspace(-Nz/2,Nz/2,Nz,q,&dq);
    dp = 2*M_PI/(Nx*dx);
    dq = 2*M_PI/(Nz*dz);
    for(int i=0; i<Nx; i++)
        p[i] *= dp;
    for(int i=0; i<Nz; i++)
        q[i] *= dq;
}

void Sim::initMatrices(){
    Phi = (std::complex<double>**)malloc(Nx*sizeof(std::complex<double>));
    for (int i=0;i<Nx;i++){
        Phi[i] = (std::complex<double>*)malloc(Nz*sizeof(std::complex<double>));
    }

    PhiMomentum = (std::complex<double>**)malloc(Nx*sizeof(std::complex<double>));
    for (int i=0;i<Nx;i++){
        PhiMomentum[i] = (std::complex<double>*)malloc(Nz*sizeof(std::complex<double>));
    }

    Phi2 = (double**)malloc(Nx*sizeof(double));
    for (int i=0;i<Nx;i++){
        Phi2[i] = (double*)malloc(Nz*sizeof(double));
    }

    PhiFuture = (std::complex<double>**)malloc(Nx*sizeof(std::complex<double>));
    for (int i=0;i<Nx;i++){
        PhiFuture[i] = (std::complex<double>*)malloc(Nz*sizeof(std::complex<double>));
    }

    cgaussian(Phi, x, z, x0, z0, a, q0, p0, Nx, Nz);
}

void Sim::planFFT(){
    in = (fftw_complex*)fftw_malloc(Nx*Nz*sizeof(fftw_complex));
    out = (fftw_complex*)fftw_malloc(Nx*Nz*sizeof(fftw_complex));
    forward = fftw_plan_dft_2d(Nx, Nz, in, out, FFTW_FORWARD, FFTW_MEASURE);
    backward = fftw_plan_dft_2d(Nx, Nz, in, out, FFTW_BACKWARD, FFTW_MEASURE);
}

void Sim::timeStep(){
    std::complex<double>* PhiTemp;
    PhiTemp = (std::complex<double>*) malloc(Nx*Nz*sizeof(std::complex<double>));
    flatten(Phi,reinterpret_cast<std::complex<double> *>(in), Nx, Nz);
    fftw_execute(forward);
    std::memcpy(in,out,sizeof(fftw_complex)*Nx*Nz);
    fftw_execute(backward);
    unflatten(reinterpret_cast<std::complex<double> *>(out), Phi, Nx, Nz);    
}
