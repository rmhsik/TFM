#include <cmath>
#include "sim.h"
#include "potential.h"

Potential::Potential(Sim *_sim){
    sim = _sim;
    omega = _sim->omega;
    l = _sim->l;
}

double Potential::envelope(double x, double z){
	if ((0.0<z) && (z<l)){
		return 1.0;
	}
	else{
		return 0.0;
	}
}

double Potential::potential(double x, double z, double t){
	return omega*envelope(x,z);

}
