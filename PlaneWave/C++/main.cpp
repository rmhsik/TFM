#include <iostream>
#include <fstream>
#include <complex>
#include <omp.h>
#include "sim.h"


#define NX      4096
#define NZ      2048
#define ZMIN    -10.0
#define ZMAX    10.0
#define XMIN    -10.0
#define XMAX    10.0
#define X0      0.0
#define Z0      0.0
#define P0      -6.0
#define Q0      -10.0
#define A       1.0

int main(){
    omp_set_num_threads(4);
    std::ofstream initWave;
    std::ofstream finalWave;
    initWave.open("InitialPhi.dat");
    finalWave.open("FinalPhi.dat");

    Sim test(XMIN,XMAX,ZMIN,
             ZMAX,X0,Z0,Q0,
             P0,A,NX,NZ);
    test.writeWavePacket(initWave);
    test.Evolution();
    test.writeWavePacket(finalWave);
    initWave.close();
    finalWave.close();
    return 0;
}