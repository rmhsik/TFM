#include <cmath>
#include <algorithm>
#include "utils.h"

void cgaussian(std::complex<double> **phi, double *x, 
                double *z, double x0, double z0, double a, 
                double q0, double p0, int Nx, int Nz){

    for(int i=0; i<Nx;i++){
        for(int j=0; j<Nz;j++){
            phi[i][j] = exp(-(pow((x[i]-x0),2)+pow(z[j]-z0,2))/(2.0*pow(a,2)))*\
                        exp(-I*p0*x[i])*exp(-I*q0*z[j]);
        }
    }
}

void waveSqr(std::complex<double> **phi, double **phi2, int Nx, int Nz){

    for(int i=0; i<Nx;i++){
        for(int j=0; j<Nz;j++){
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
    for(int i=0;i<Nx;i++){
        for(int j=0; j<Nz; j++){
            out[i*Nz+j] = in[i][j];
        }
    }
}

void unflatten(std::complex<double> *in, std::complex<double> **out, int Nx, int Nz){
    for(int i=0;i<Nx;i++){
        for(int j=0; j<Nz; j++){
            out[i][j] = in[i*Nz+j];
        }
    }
}

void normalize(std::complex<double> **in, int Nx, int Nz){
    for(int i=0; i < Nx; i++){
        for(int j = 0; j<Nz; j++){
            in[i][j] /= (Nx*Nz);
        }
    }
}

void freqshift(double *in, double *out, int N){
    std::copy(in+N/2,in+N,out);
    std::copy(in, in+N/2, out + N/2);
}