#include <iostream>
#include <cstring>
#include "evolution.h"
#include "potential.h"
#include "utils.h"
#include "sim.h"


Evolution::Evolution(Sim *_sim){
    sim = _sim;
    forward = sim->forward;
    backward = sim->backward;
    fftwin = sim->in;
    fftwout = sim->out;

    temp1 = (std::complex<double> *)malloc(sim->Nx*sim->Nz*sizeof(std::complex<double>));
    temp2 = (std::complex<double> *)malloc(sim->Nx*sim->Nz*sizeof(std::complex<double>));
}

void Evolution::evMomentum(std::complex<double> *in,std::complex<double> *out,  double t){
    int i, j;
    #pragma omp parallel for private(i, j)
    for(i = 0; i<sim->Nz; i++){
        for (j = 0; j<sim->Nx; j++){
            out[i*sim->Nx + j] = in[i*sim->Nx+j]*exp(-I*pow(sim->pshift[j],2)/(2*sim->m)*sim->dt)*\
                            exp(-I*pow(sim->qshift[i],2)/(2*sim->m)*sim->dt);
        }
    }
}

void Evolution::evSpace(std::complex<double> *in1, std::complex<double> *in2, std::complex<double> *out, double dt){
	int i,j;

	#pragma omp parallel for private(i, j)
	for (i = 0; i < sim->Nz; i++){
		for(j = 0; j < sim->Nx; j++){
	    		out[i*sim->Nx + j] = cos(sim->potOperator->potential(sim->x[i],sim->z[i],0.0)*dt)*in1[i*sim->Nx + j]
                                  -I*sin(sim->potOperator->potential(sim->x[i],sim->z[i],0.0)*dt)*in2[i*sim->Nx + j];
	    }
	}

}

/*
void Evolution::evSpace(std::complex<double> *in, std::complex<double> *out, double dt){
	double V0 = 0.0;
	int i,j;

	#pragma omp parallel for private(i, j)
	for (i = 0; i < sim->Nz; i++){
		for(j = 0; j < sim->Nx; j++){
			if((1250<i) && (i<1300))
	    		out[i*sim->Nx + j] = in[i*sim->Nx + j]*exp(-I*V0*dt);
    
			out[i*sim->Nx + j] = in[i*sim->Nx + j];
	    }
	}

}
*/

void Evolution::timeStepG(double t){
    flatten(sim->PhiG, temp1, sim->Nx, sim->Nz);
    flatten(sim->PhiE, temp2, sim->Nx, sim->Nz);

    evSpace(temp1, temp2, reinterpret_cast<std::complex<double> *>(fftwin), sim->dt);
    fftw_execute(forward);

    evMomentum(reinterpret_cast<std::complex<double> *>(fftwout),reinterpret_cast<std::complex<double> *>(fftwout),sim->dt);
    std::memcpy(fftwin,fftwout,sizeof(fftw_complex)*sim->Nx*sim->Nz);
    fftw_execute(backward);
    
    unflatten(reinterpret_cast<std::complex<double> *>(fftwout), sim->PhiFutureG, sim->Nx, sim->Nz);
    //normalize(sim->PhiG,sim->Nx,sim->Nz);
}

void Evolution::timeStepE(double t){
    flatten(sim->PhiG, temp1, sim->Nx, sim->Nz);
    flatten(sim->PhiE, temp2, sim->Nx, sim->Nz);

    evSpace(temp2, temp1, reinterpret_cast<std::complex<double> *>(fftwin), sim->dt);
    fftw_execute(forward);

    evMomentum(reinterpret_cast<std::complex<double> *>(fftwout),reinterpret_cast<std::complex<double> *>(fftwout),sim->dt);
    std::memcpy(fftwin,fftwout,sizeof(fftw_complex)*sim->Nx*sim->Nz);
    fftw_execute(backward);
    
    unflatten(reinterpret_cast<std::complex<double> *>(fftwout), sim->PhiFutureE, sim->Nx, sim->Nz);
    //normalize(sim->PhiE,sim->Nx,sim->Nz);
}

///////////////////////////////////////////////////////////

void Evolution::timeStep(double t){
    timeStepG(t);
    timeStepE(t);

    normalize(sim->PhiFutureG,sim->PhiG,sim->Nx, sim-> Nz);
    normalize(sim->PhiFutureE,sim->PhiE,sim->Nx, sim-> Nz);

}