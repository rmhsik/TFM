#include <iostream>
#include <complex>
#include <omp.h>
#include "sim.h"
#include "config.h"

int main(){
    omp_set_num_threads(4);
    //std::ofstream initWave;
    //std::ofstream finalWave;
    //initWave.open("InitialPhi.dat");
    //finalWave.open("FinalPhi.dat");
    Config conf;

    //Sim test(conf.XMIN,conf.XMAX,conf.ZMIN,
    //         conf.ZMAX,conf.TMIN,conf.TMAX,conf.X0,conf.Z0,conf.Q0,
    //         conf.P0,conf.A,conf.NX,conf.NZ,conf.NT,conf.NSAMPLE);
    Sim test(&conf);
    //test.writeWavePacket(initWave);
    test.Run();
    //test.writeWavePacket(finalWave);
    //initWave.close();
    //finalWave.close();
    return 0;
}
