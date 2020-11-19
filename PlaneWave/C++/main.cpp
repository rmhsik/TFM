#include <iostream>
#include <complex>
#include <omp.h>
#include "sim.h"


#define NX      4096
#define NZ      2048
#define NT      1000
#define ZMIN    -10.0
#define ZMAX    10.0
#define XMIN    -10.0
#define XMAX    10.0
#define TMAX    0.7
#define TMIN    0.0
#define X0      2.5
#define Z0      2.5
#define P0      -6.0
#define Q0      -10.0
#define A       0.5
#define NSAMPLE 10

int main(){
    omp_set_num_threads(4);
    //std::ofstream initWave;
    //std::ofstream finalWave;
    //initWave.open("InitialPhi.dat");
    //finalWave.open("FinalPhi.dat");

    Sim test(XMIN,XMAX,ZMIN,
             ZMAX,TMIN,TMAX,X0,Z0,Q0,
             P0,A,NX,NZ,NT,NSAMPLE);
    //test.writeWavePacket(initWave);
    test.Run();
    //test.writeWavePacket(finalWave);
    //initWave.close();
    //finalWave.close();
    return 0;
}
