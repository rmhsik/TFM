#include "sim.h"
#include "potential.h"

Potential::Potential(Sim *_sim){
    sim = _sim;
    omega = _sim->omega;
}