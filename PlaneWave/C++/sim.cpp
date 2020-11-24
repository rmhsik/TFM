#include <iostream>
#include <fstream>
#include <cstring>
#include <chrono>
#include <thread>
#include "sim.h"
#include "utils.h"
#include "evolution.h"
#include "config.h"
#include "potential.h"

Sim::Sim(Config *_conf){
    
    Nx = _conf->NX;
    Nz = _conf->NZ;
    Nt = _conf->NT;
    xmin = _conf->XMIN;
    xmax = _conf->XMAX;
    zmin = _conf->ZMIN;
    zmax = _conf->ZMAX;
    tmin = _conf->TMIN;
    tmax = _conf->TMAX;
    x0 = _conf->X0;
    z0 = _conf->Z0;
    q0 = _conf->Q0;
    p0 = _conf->P0;
    a = _conf->A;
    Nsample = _conf->NSAMPLE;
    m = _conf->M;
    l = _conf->L;
    omega = _conf->OMEGA;

    auto t1 = std::chrono::high_resolution_clock::now();
    initSpaceTime();
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

    evOperator = new Evolution(this);
    conf = _conf;
	potOperator = new Potential(this);
}

void Sim::writeWavePacket(int j){
    const char nameG[5] = "PhiG";
    const char nameE[5] = "PhiE";
    waveSqr(PhiG, Phi2G, Nx, Nz);
    waveSqr(PhiE, Phi2E, Nx, Nz);
    std::thread writeG(&Sim::write2File,this,PhiG,Phi2G,nameG,j);
    std::thread writeE(&Sim::write2File,this,PhiE,Phi2E,nameE,j);
    std::thread writeT(&Sim::write2FileT,this,j);
    writeG.join();
    writeE.join();
    writeT.join();
}

void Sim::Benchmark(){
    auto t1 = std::chrono::high_resolution_clock::now();
    evOperator->timeStep(0.0);
    auto t2 = std::chrono::high_resolution_clock::now();    
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();
    std::cout<<"Estimated time for simulation: " << Nt*duration<<std::endl;
}

void Sim::Run(){
    const int DeltaN = Nt/Nsample;
    int j = 1;
    for(int i=0; i<Nt; i++){
        evOperator->timeStep(t[i]);
        //std::cout<<i<<std::endl;
        if(i%DeltaN == 0){
            std::cout<<"Saving "<<j<<" of "<< Nsample<<std::endl;
            writeWavePacket(i);
            j++;
        }
    }
}

///////////////////////////////////////////////////////////////////

void Sim::initSpaceTime(){
    x = (double*)malloc(Nx*sizeof(double));
    z = (double*)malloc(Nz*sizeof(double));
    t = (double*)malloc(Nt*sizeof(double));
    linspace(xmin, xmax, Nx, x, &dx);
    linspace(zmin, zmax, Nz, z, &dz);
    linspace(tmin, tmax, Nt, t, &dt);
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
    int i;
    PhiG = (std::complex<double>**)malloc(Nz*sizeof(std::complex<double>));
    #pragma omp parallel for private(i)
    for (i=0;i<Nz;i++){
        PhiG[i] = (std::complex<double>*)malloc(Nx*sizeof(std::complex<double>));
    }

    PhiMomentumG = (std::complex<double>**)malloc(Nz*sizeof(std::complex<double>));
    #pragma omp parallel for private(i)
    for (i=0;i<Nz;i++){
        PhiMomentumG[i] = (std::complex<double>*)malloc(Nx*sizeof(std::complex<double>));
    }

    Phi2G = (double**)malloc(Nz*sizeof(double));
    #pragma omp parallel for private(i)
    for (i=0;i<Nz;i++){
        Phi2G[i] = (double*)malloc(Nx*sizeof(double));
    }

    PhiFutureG = (std::complex<double>**)malloc(Nz*sizeof(std::complex<double>));
    #pragma omp parallel for private(i)
    for (i=0;i<Nz;i++){
        PhiFutureG[i] = (std::complex<double>*)malloc(Nx*sizeof(std::complex<double>));
    }

    PhiE = (std::complex<double>**)malloc(Nz*sizeof(std::complex<double>));
    #pragma omp parallel for private(i)
    for (i=0;i<Nz;i++){
        PhiE[i] = (std::complex<double>*)malloc(Nx*sizeof(std::complex<double>));
    }

    PhiMomentumE = (std::complex<double>**)malloc(Nz*sizeof(std::complex<double>));
    #pragma omp parallel for private(i)
    for (i=0;i<Nz;i++){
        PhiMomentumE[i] = (std::complex<double>*)malloc(Nx*sizeof(std::complex<double>));
    }

    Phi2E = (double**)malloc(Nz*sizeof(double));
    #pragma omp parallel for private(i)
    for (i=0;i<Nz;i++){
        Phi2E[i] = (double*)malloc(Nx*sizeof(double));
    }

    PhiFutureE = (std::complex<double>**)malloc(Nz*sizeof(std::complex<double>));
    #pragma omp parallel for private(i)
    for (i=0;i<Nz;i++){
        PhiFutureE[i] = (std::complex<double>*)malloc(Nx*sizeof(std::complex<double>));
    }
    
    cgaussian(PhiG, x, z, x0, z0, a, q0, p0, Nx, Nz);
    cgaussian(PhiE, x, z, x0, z0, a, 0.0, 0.0, Nx, Nz);
}

void Sim::planFFT(){
    fftw_init_threads();
    fftw_plan_with_nthreads(4);
    in = (fftw_complex*)fftw_malloc(Nx*Nz*sizeof(fftw_complex));
    out = (fftw_complex*)fftw_malloc(Nx*Nz*sizeof(fftw_complex));
    forward = fftw_plan_dft_2d(Nz, Nx, in, out, FFTW_FORWARD, FFTW_MEASURE);
    backward = fftw_plan_dft_2d(Nz, Nx, in, out, FFTW_BACKWARD, FFTW_MEASURE);
}

void Sim::write2File(std::complex<double> **Phi,double **Phi2,
                    const char *name, int j){
    char path[60];
    char step[5];
    strcpy(path,"Data/");
    strcat(path,name);
    strcat(path,"/");
    strcat(path,name);
    sprintf(step,"%d",j);
    strcat(path,step);
    strcat(path,".dat");
    std::ofstream file;
    file.open(path);

    //waveSqr(Phi, Phi2, Nx, Nz);
    for (int i=0;i<Nz; i++){
        for (int j=0; j<Nx; j++){
            file << Phi2[i][j]<<" ";
        }
        file<<"\n";
    }
    file.close();
}

void Sim::write2FileT(int j){
    char path[60];
    char step[5];
    strcpy(path,"Data/PhiT/PhiT");
    sprintf(step,"%d",j);
    strcat(path,step);
    strcat(path,".dat");
    std::ofstream file;
    file.open(path);

    for (int i=0;i<Nz; i++){
        for (int j=0; j<Nx; j++){
            file << 0.5*(Phi2G[i][j]+Phi2E[i][j])<<" ";
        }
        file<<"\n";
    }
    file.close();
}
