#include <cmath>
#include <algorithm>
#include "utils.h"

void cgaussian(std::complex<double> **phi, double *x, 
                double *z, double x0, double z0, double a, 
                double q0, double p0, int Nx, int Nz){
    int i,j;
    #pragma omp parallel for private(i, j)
    for(i=0; i<Nz;i++){
        for(j=0; j<Nx;j++){
            phi[i][j] = exp(-(pow((x[j]-x0),2)+pow(z[i]-z0,2))/(2.0*pow(a,2)))*\
                        exp(I*p0*x[j])*exp(I*q0*z[i]);
        }
    }
}

void zeros(std::complex<double> **phi, int Nx, int Nz){
    int i, j;
    #pragma omp parallel for private(i, j)
    for(i=0; i<Nz;i++){
        for(j=0; j<Nx;j++){
            phi[i][j] = 0.0;
        }
    }
}

void waveSqr(std::complex<double> **phi, double **phi2, int Nx, int Nz){

    for(int i=0; i<Nz;i++){
        for(int j=0; j<Nx;j++){
            phi2[i][j] = pow(abs(phi[i][j]),2);
        }
    }
}


void linspace(double a, double b, int N, double *out, double *dz){
    *dz = fabs(a-b)/(double)N;
    for (int i = 0; i < N; i++){
        out[i] = a + *dz*i; 
    }
}

void flatten(std::complex<double> **in, std::complex<double> *out, int Nx, int Nz){
    for(int i=0;i<Nz;i++){
        for(int j=0; j<Nx; j++){
            out[i*Nx+j] = in[i][j];
        }
    }
}

void unflatten(std::complex<double> *in, std::complex<double> **out, int Nx, int Nz){
    for(int i=0;i<Nz;i++){
        for(int j=0; j<Nx; j++){
            out[i][j] = in[i*Nx+j];
        }
    }
}

void normalize(std::complex<double> **in, std::complex<double> **out, int Nx, int Nz){
    for(int i=0; i < Nz; i++){
        for(int j = 0; j<Nx; j++){
            out[i][j] /= (Nx*Nz);
        }
    }
}

void freqshift(double *in, double *out, int N){
    std::copy(in+N/2,in+N,out);
    std::copy(in, in+N/2, out + N/2);
}