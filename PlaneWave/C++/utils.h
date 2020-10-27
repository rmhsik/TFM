#ifndef UTILS_H
#define UILS_H
#include <complex>
#include <fftw3.h>

const static std::complex<double> I(std::complex<double>(0,1));

void cgaussian(std::complex<double> **phi, double *x, 
                double *z, double x0, double z0, double a, 
                double q0, double p0, int Nx, int Nz);

void waveSqr(std::complex<double> **phi,double **phi2, int Nx, int Nz);

void linspace(double a, double b, int N, double *out, double *dz);

void flatten(std::complex<double> **in, std::complex<double> *out, int Nx, int Nz);

void unflatten(std::complex<double> *in, std::complex<double> **out, int Nx, int Nz);

#endif
