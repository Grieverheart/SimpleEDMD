#include <iostream>
#include "include/Time.h"
#include "include/EventManager.h"
#include "include/Simulation.h"


int main(int argc, char *argv[]){

    Simulation sim;

    sim.readConfig(argv[1]);

    sim.init();

    sim.run();

    return 0;
}
