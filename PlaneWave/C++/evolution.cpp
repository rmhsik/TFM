#include <iostream>
#include <fstream>
#include <cstring>
#include <complex>
#include <fftw3.h>
#include <cmath>
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
    omega = _sim->omega;
    l = _sim->l;

    temp1 = (std::complex<double> *)malloc(sim->Nx*sim->Nz*sizeof(std::complex<double>));
    temp2 = (std::complex<double> *)malloc(sim->Nx*sim->Nz*sizeof(std::complex<double>));

    std::ofstream file;
    file.open("Data/Potential.dat");
    std::cout<<sin(potential(4.0,1.5, 0.0)*sim->dt)<<std::endl;
    for (int i = 0; i < sim->Nz; i++){
		for(int j = 0; j < sim->Nx; j++){
            //if (sim->potOperator->potential(sim->x[i],sim->z[j],0.0) > 0.0)
            file<<sin(potential(sim->x[j],sim->z[i], 0.0)*sim->dt)<<" ";
        }
        file<<std::endl;
    }
    file.close();
}

double Evolution::envelope(double x, double z){
	//if ((0.0<z) && (z<l)){
		return omega;
	//}
	//else{
	//	return 0.0;
	//}
}

double Evolution::potential(double x, double z, double t){
	
	return envelope(x,z);

}


void Evolution::evMomentumG(std::complex<double> *in,std::complex<double> *out,  double t){
    int i, j;
    #pragma omp parallel for private(i, j)
    for(i = 0; i<sim->Nz; i++){
        for (j = 0; j<sim->Nx; j++){
            out[i*sim->Nx + j] = in[i*sim->Nx+j]*exp(-I*pow(sim->pshift[j],2)/(2*sim->m)*sim->dt)*\
                            exp(-I*pow(sim->qshift[i],2)/(2*sim->m)*sim->dt);
        }
    }
}

void Evolution::evMomentumE(std::complex<double> *in,std::complex<double> *out,  double t){
    int i, j;
    #pragma omp parallel for private(i, j)
    for(i = 0; i<sim->Nz; i++){
        for (j = 0; j<sim->Nx; j++){
            out[i*sim->Nx + j] = in[i*sim->Nx+j]*exp(-I*pow(sim->pshift[j]+sim->k,2)/(2*sim->m)*sim->dt)*\
                            exp(-I*pow(sim->qshift[i],2)/(2*sim->m)*sim->dt);
        }
    }
}

void Evolution::evSpaceG(std::complex<double> *in1, std::complex<double> *in2, std::complex<double> *out, double dt){
	int i,j;

	#pragma omp parallel for private(i, j)
	for (i = 0; i < sim->Nz; i++){
		for(j = 0; j < sim->Nx; j++){
                //std::cout<<-I*sin(sim->potOperator->potential(sim->x[i],sim->z[i],0.0)*dt)*in2[i*sim->Nx + j]<<std::endl;
	    		out[i*sim->Nx + j] = cos(potential(sim->x[j],sim->z[i],0.0)*dt)*in1[i*sim->Nx + j]
                                  -I*sin(potential(sim->x[j],sim->z[i],0.0)*dt)*in2[i*sim->Nx + j];
	    }
	}

}

void Evolution::evSpaceE(std::complex<double> *in1, std::complex<double> *in2, std::complex<double> *out, double dt){
	int i,j;

	#pragma omp parallel for private(i, j)
	for (i = 0; i < sim->Nz; i++){
		for(j = 0; j < sim->Nx; j++){
                //std::cout<<-I*sin(sim->potOperator->potential(sim->x[i],sim->z[i],0.0)*dt)*in2[i*sim->Nx + j]<<std::endl;
	    		out[i*sim->Nx + j] = cos(potential(sim->x[j],sim->z[i],0.0)*dt)*in2[i*sim->Nx + j]
                                  -I*sin(potential(sim->x[j],sim->z[i],0.0)*dt)*in1[i*sim->Nx + j];

                //if (pow(abs(out[i*sim->Nx + j]),2)>1E-1)
                //    std::cout<<pow(abs(out[i*sim->Nx + j]),2)<<std::endl;
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

    evSpaceG(temp1, temp2, reinterpret_cast<std::complex<double> *>(fftwin), sim->dt);
    //fftw_execute(forward);

    //evMomentumG(reinterpret_cast<std::complex <double> *>(fftwout),reinterpret_cast<std::complex<double> *>(fftwin),sim->dt);
    //std::memcpy(fftwin,fftwout,sizeof(fftw_complex)*sim->Nx*sim->Nz);
    //fftw_execute(backward);
    
    unflatten(reinterpret_cast<std::complex<double> *>(fftwout), sim->PhiFutureG, sim->Nx, sim->Nz);
    //normalize(sim->PhiFutureG,sim->PhiFutureG,sim->Nx, sim-> Nz);
}

void Evolution::timeStepE(double t){
    flatten(sim->PhiG, temp1, sim->Nx, sim->Nz);
    flatten(sim->PhiE, temp2, sim->Nx, sim->Nz);

    evSpaceE(temp1, temp2, reinterpret_cast<std::complex<double> *>(fftwin), sim->dt);
    //fftw_execute(forward);

    //evMomentumE(reinterpret_cast<std::complex<double> *>(fftwout),reinterpret_cast<std::complex<double> *>(fftwin),sim->dt);
   // std::memcpy(fftwin,fftwout,sizeof(fftw_complex)*sim->Nx*sim->Nz);
    //fftw_execute(backward);
    
    unflatten(reinterpret_cast<std::complex<double> *>(fftwout), sim->PhiFutureE, sim->Nx, sim->Nz);
    //normalize(sim->PhiFutureE,sim->PhiFutureE,sim->Nx,sim->Nz);
}

///////////////////////////////////////////////////////////

void Evolution::timeStep(double t){
    timeStepE(t);
    timeStepG(t);

    normalize(sim->PhiFutureG,sim->PhiG, 1, 1);
    normalize(sim->PhiFutureE,sim->PhiE, 1, 1);

}