#ifndef __SIMULATION_H
#define __SIMULATION_H

#include <vector>
#include <random>
#include "BoundaryCondition.h"
#include "Vec.h"
#include "EventManager.h"
#include "CellList.h"
#include "Particle.h"

//NOTE:
//1. The Boundary Condition class will most probably need to be handed
//   as a template parameter for performance reasons.
//2. Make a Collision Response class for different collision responses.
//   We might have to pass this as a template parameter too.
//3. We might also need to switch to -B/2 -> +B/2 coordinates so that
//   we can easily implement spherical boundary conditions.
//4. We should have two ways to create a Simulation. One with a
//   fromFile method and one by handing the needed information.

class Simulation{
public:
    Simulation(const CubicPBC& pbc, std::vector<Particle>&& particles):
        nSpheres_(particles.size()), time_(0.0),
        pbc_(pbc), particles_(particles)
    {
        mtGen_.seed(0);
    }

    void run(void);
    bool init(void);
    void saveConfig(const char* filename);
private:
    ParticleEvent getCollisionEvent(int pA, int pB)const;
    ParticleEvent getCellCrossEvent(int pid)const;
    void runCollisionEvent(const ParticleEvent& event);
    void runCellCrossEvent(const ParticleEvent& event);
    void updateParticle(int pid);

    int    nSpheres_;
    double time_;

    CubicPBC pbc_;

    std::vector<int>      nCollisions_;
    std::vector<Particle> particles_;

    EventManager eventManager_;
    CellList     cll_;

    std::mt19937 mtGen_;
};

#endif
