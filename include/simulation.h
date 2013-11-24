#ifndef __SIMULATION_H
#define __SIMULATION_H

#include "Vec.h"

//Start by creating a simple simulation class so we can start testing.
//After sufficient testing, we can move to the entity-component system
//and check how performance is affected. If the performance penalty is
//not too significant, it will beneficial to keep the EC system.

class Simulation{
public:
    void run(void);
    void addSphere(Vec3d position, double radius);
private:
    int nSpheres;
    std::vector<double> radii;
    std::vector<Vec3d>  positions;
};

#endif
