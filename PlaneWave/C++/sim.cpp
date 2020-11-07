#include <iostream>
#include <fstream>
#include <cstring>
#include <chrono>
#include "sim.h"
#include "utils.h"

Sim::Sim(double _xmin, double _xmax, double _zmin,
         double _zmax, double _tmin, double _tmax,
         double _x0, double _z0,
         double _q0, double _p0, 
         double _a, int _Nx, int _Nz , int _Nt, int _Nsample){
    
    Nx = _Nx;
    Nz = _Nz;
    Nt = _Nt;
    xmin = _xmin;
    xmax = _xmax;
    zmin = _zmin;
    zmax = _zmax;
    tmin = _tmin;
    tmax = _tmax;
    x0 = _x0;
    z0 = _z0;
    q0 = _q0;
    p0 = _p0;
    a = _a;
    Nsample = _Nsample;

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
}

void Sim::writeWavePacket(int j){
    char path[60];
    char step[5];
    strcpy(path,"Data/Phi");
    sprintf(step,"%d",j);
    strcat(path,step);
    strcat(path,".dat");
    std::ofstream file;
    file.open(path);
    waveSqr(Phi, Phi2, Nx, Nz);
    for (int i=0;i<Nz; i++){
        for (int j=0; j<Nx; j++){
            file << Phi2[i][j]<<" ";
        }
        file<<"\n";
    }
    file.close();
}

void Sim::Benchmark(){
    auto t1 = std::chrono::high_resolution_clock::now();
    timeStep(0.0);
    auto t2 = std::chrono::high_resolution_clock::now();    
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();
    std::cout<<"Estimated time for simulation: " << Nt*duration<<std::endl;
}

void Sim::Evolution(){
    const int DeltaN = Nt/Nsample;
    int j = 1;
    for(int i=0; i<Nt; i++){
        timeStep(t[i]);
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

void Sim::evMomentum(std::complex<double> *in, double dt){
    double m = 1;
    int i, j;
    #pragma omp parallel for private(i, j)
    for(i = 0; i<Nz; i++){
        for (j = 0; j<Nx; j++){
            in[i*Nx + j] = in[i*Nx+j]*exp(-I*pow(pshift[j],2)/(2*m)*dt)*\
                            exp(-I*pow(qshift[i],2)/(2*m)*dt);
        }
    }

}

void Sim::evSpace(std::complex<double> *in, double dt){
	double V0 = 40.0;
	int i,j;

	#pragma omp parallel for private(i, j)
	for (i = 0; i < Nz; i++){
		for(j = 0; j < Nx; j++){
			if((1250<i) && (i<1270))
				in[i*Nx + j] = in[i*Nx + j]*exp(-I*V0*dt);

			in[i*Nx + j] = in[i*Nx + j];
		}
	}

}


void Sim::timeStep(double t){
    //std::cout<<dt<<std::endl;
    flatten(Phi,reinterpret_cast<std::complex<double> *>(in), Nx, Nz);
    evSpace(reinterpret_cast<std::complex<double> *>(in),dt/2.0);
    fftw_execute(forward);
    evMomentum(reinterpret_cast<std::complex<double> *>(out),dt);
    std::memcpy(in,out,sizeof(fftw_complex)*Nx*Nz);
    fftw_execute(backward);
	evSpace(reinterpret_cast<std::complex<double> *>(out),dt/2.0);
    unflatten(reinterpret_cast<std::complex<double> *>(out), Phi, Nx, Nz);
    normalize(Phi,Nx,Nz);
}

