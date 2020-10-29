#include <iostream>
#include <cstring>
#include <chrono>
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

    auto t1 = std::chrono::high_resolution_clock::now();
    initSpace();
    initMomentum();
    initMatrices();
    planFFT();
    auto t2 = std::chrono::high_resolution_clock::now();    
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();
    std::cout<<"Time used for intiaizing: " << duration<<std::endl;

    qshift = (double*)malloc(Nz*sizeof(double));
    pshift = (double*)malloc(Nx*sizeof(double));
    freqshift(p, pshift, Nx);
    freqshift(q, qshift, Nz);
}

void Sim::writeWavePacket(std::ofstream& file){
    waveSqr(Phi, Phi2, Nx, Nz);
    for (int i=0;i<Nz; i++){
        for (int j=0; j<Nx; j++){
            file << Phi2[i][j]<<" ";
        }
        file<<"\n";
    }
}

void Sim::Evolution(){
    auto t1 = std::chrono::high_resolution_clock::now();
    timeStep();
    auto t2 = std::chrono::high_resolution_clock::now();    
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();
    std::cout<<"Time used for one timestep: " << duration<<std::endl;
}

///////////////////////////////////////////////////////////////////

void Sim::initSpace(){
    x = (double*)malloc(Nx*sizeof(double));
    z = (double*)malloc(Nz*sizeof(double));
    linspace(xmin,xmax,Nx,x,&dx);
    linspace(zmin,zmax,Nz,z,&dz);
}

void Sim::initMomentum(){
    q = (double*)malloc(Nz*sizeof(double));
    p = (double*)malloc(Nx*sizeof(double));

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
    Phi = (std::complex<double>**)malloc(Nz*sizeof(std::complex<double>));
    #pragma omp parallel for
    for (int i=0;i<Nz;i++){
        Phi[i] = (std::complex<double>*)malloc(Nx*sizeof(std::complex<double>));
    }

    PhiMomentum = (std::complex<double>**)malloc(Nz*sizeof(std::complex<double>));
    #pragma omp parallel for
    for (int i=0;i<Nz;i++){
        PhiMomentum[i] = (std::complex<double>*)malloc(Nx*sizeof(std::complex<double>));
    }

    Phi2 = (double**)malloc(Nz*sizeof(double));
    #pragma omp parallel for
    for (int i=0;i<Nz;i++){
        Phi2[i] = (double*)malloc(Nx*sizeof(double));
    }

    PhiFuture = (std::complex<double>**)malloc(Nz*sizeof(std::complex<double>));
    #pragma omp parallel for
    for (int i=0;i<Nz;i++){
        PhiFuture[i] = (std::complex<double>*)malloc(Nx*sizeof(std::complex<double>));
    }

    cgaussian(Phi, x, z, x0, z0, a, q0, p0, Nx, Nz);
}

void Sim::planFFT(){
    fftw_init_threads();
    fftw_plan_with_nthreads(4);
    in = (fftw_complex*)fftw_malloc(Nx*Nz*sizeof(fftw_complex));
    out = (fftw_complex*)fftw_malloc(Nx*Nz*sizeof(fftw_complex));
    forward = fftw_plan_dft_2d(Nz, Nx, in, out, FFTW_FORWARD, FFTW_MEASURE);
    backward = fftw_plan_dft_2d(Nz, Nx, in, out, FFTW_BACKWARD, FFTW_MEASURE);
}

void Sim::evMomentum(std::complex<double> *in, double t){
    double m = 1;
    int i, j;
    #pragma omp parallel for private(i, j)
    for(i = 0; i<Nz; i++){
        for (j = 0; j<Nx; j++){
            in[i*Nx + j] = in[i*Nx+j]*exp(-I*pow(pshift[j],2)/(2*m)*t)*\
                            exp(-I*pow(qshift[i],2)/(2*m)*t);
        }
    }

}

void Sim::timeStep(){
    double t = 0.7;
    flatten(Phi,reinterpret_cast<std::complex<double> *>(in), Nx, Nz);
    fftw_execute(forward);
    evMomentum(reinterpret_cast<std::complex<double> *>(out),t);
    std::memcpy(in,out,sizeof(fftw_complex)*Nx*Nz);
    fftw_execute(backward);
    unflatten(reinterpret_cast<std::complex<double> *>(out), Phi, Nx, Nz);
    normalize(Phi,Nx,Nz);
}

