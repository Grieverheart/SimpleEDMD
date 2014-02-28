#ifndef __SIMULATION_H
#define __SIMULATION_H

#include <vector>
#include <random>
#include "Vec.h"
#include "EventManager.h"
#include "Particle.h"

//NOTE:
//1. Make Boundary Conditions a class so that we can change it easily.
//2. Make a Collision Response class for different collision responses.
//3. Will probably introduce the notion of properties for quantities such
//   as mass.

//NOTE: Consider using smart pointers for passing events around
class Simulation{
public:
    Simulation(void):
        nSpheres_(0), time_(0.0)
    {
        mtGen_.seed(0);
    }

    void run(void);
    bool init(void);
    void readConfig(const char* filename);
    void saveConfig(const char* filename);
private:
    CollisionEvent* getCollisionEvent(size_t pA, size_t pB)const;
    void runCollisionEvent(const CollisionEvent& event);
    void updateParticle(size_t pID);

    Vec3d applyPeriodicBC(const Vec3d& vec)const;

    size_t nSpheres_;
    double boxSize_;
    double time_;

    std::vector<size_t>   nCollisions_;
    std::vector<Particle> particles_;
    std::vector<double>   radii_;

    EventManager eventManager_;

    std::mt19937 mtGen_;
};

#endif
