#include <iostream>
#include <cstring>
#include "evolution.h"
#include "utils.h"


Evolution::Evolution(Sim *_sim){
    sim = _sim;
    forward = sim->forward;
    backward = sim->backward;
    fftwin = sim->in;
    fftwout = sim->out;
}

void Evolution::evMomentum(std::complex<double> *in, double t){
    int i, j;
    #pragma omp parallel for private(i, j)
    for(i = 0; i<sim->Nz; i++){
        for (j = 0; j<sim->Nx; j++){
            in[i*sim->Nx + j] = in[i*sim->Nx+j]*exp(-I*pow(sim->pshift[j],2)/(2*sim->m)*sim->dt)*\
                            exp(-I*pow(sim->qshift[i],2)/(2*sim->m)*sim->dt);
        }
    }
}

void Evolution::evSpace(std::complex<double> *in, double dt){
	double V0 = 0.0;
	int i,j;

	#pragma omp parallel for private(i, j)
	for (i = 0; i < sim->Nz; i++){
		for(j = 0; j < sim->Nx; j++){
			if((1250<i) && (i<1300))
	    		in[i*sim->Nx + j] = in[i*sim->Nx + j]*exp(-I*V0*dt);
    
			in[i*sim->Nx + j] = in[i*sim->Nx + j];
	    }
	}

}
void Evolution::timeStepG(double t){
    //std::cout<<dt<<std::endl;
    flatten(sim->PhiG,reinterpret_cast<std::complex<double> *>(fftwin), sim->Nx, sim->Nz);
    evSpace(reinterpret_cast<std::complex<double> *>(fftwin),sim->dt);
    fftw_execute(forward);
    evMomentum(reinterpret_cast<std::complex<double> *>(fftwout),sim->dt);
    std::memcpy(fftwin,fftwout,sizeof(fftw_complex)*sim->Nx*sim->Nz);
    fftw_execute(backward);
	//evSpace(reinterpret_cast<std::complex<double> *>(fftwout),dt/2.0);
    unflatten(reinterpret_cast<std::complex<double> *>(fftwout), sim->PhiG, sim->Nx, sim->Nz);
    normalize(sim->PhiG,sim->Nx,sim->Nz);
}

void Evolution::timeStepE(double t){
    //std::cout<<dt<<std::endl;
    flatten(sim->PhiE,reinterpret_cast<std::complex<double> *>(fftwin), sim->Nx, sim->Nz);
    evSpace(reinterpret_cast<std::complex<double> *>(fftwin),sim->dt);
    fftw_execute(forward);
    evMomentum(reinterpret_cast<std::complex<double> *>(fftwout),sim->dt);
    std::memcpy(fftwin,fftwout,sizeof(fftw_complex)*sim->Nx*sim->Nz);
    fftw_execute(backward);
	//evSpace(reinterpret_cast<std::complex<double> *>(fftwout),dt/2.0);
    unflatten(reinterpret_cast<std::complex<double> *>(fftwout), sim->PhiE, sim->Nx, sim->Nz);
    normalize(sim->PhiE,sim->Nx,sim->Nz);
}

///////////////////////////////////////////////////////////

void Evolution::timeStep(double t){
    timeStepG(t);
    timeStepE(t);
}