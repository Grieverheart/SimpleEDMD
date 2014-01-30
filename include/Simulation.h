#ifndef __SIMULATION_H
#define __SIMULATION_H

#include <vector>
#include "Vec.h"

//Start by creating a simple simulation class so we can start testing.
//After sufficient testing, we can move to the entity-component system
//and check how performance is affected. If the performance penalty is
//not too significant, it will beneficial to keep the EC system.

//NOTE: The plan is, for a particle i, to only schedule the next collision
//and eventually, when the cell list is implemented, the next cell transfer.
//check 'Algorithms for Particle-Field Simulations with Collisions'

class Simulation{
public:
    Simulation(void){};
    void run(void);
    void addSphere(Vec3d position, double radius);
    void readConfig(const char* filename);
private:
    int nSpheres_;
    double boxSize_;
    std::vector<double> radii_;
    std::vector<Vec3d>  positions_;
};

#endif
