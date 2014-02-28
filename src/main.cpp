#include <iostream>
#include "include/Simulation.h"


int main(int argc, char *argv[]){

    Simulation sim;

    sim.readConfig(argv[1]);

    sim.init();

    sim.run();

    return 0;
}
