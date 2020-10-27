#include <iostream>
#include <fstream>
#include <complex>
#include "sim.h"

#define NX      4096
#define NZ      4096
#define ZMIN    -10.0
#define ZMAX    10.0
#define XMIN    -10.0
#define XMAX    10.0
#define X0      0.0
#define Z0      0.0
#define P0      6.0
#define Q0      -10.0
#define A       1.0

int main(){
    std::ofstream initWave;
    std::ofstream finalWave;
    initWave.open("InitialPhi.dat");
    finalWave.open("FinalPhi.dat");

    Sim test(XMIN,XMAX,ZMIN,ZMAX,X0,Z0,Q0,P0,A,NZ,NZ);
    test.writeWavePacket(initWave);
    test.Evolution();
    test.writeWavePacket(finalWave);
    initWave.close();
    finalWave.close();
    return 0;
}