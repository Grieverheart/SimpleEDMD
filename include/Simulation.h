#ifndef __SIMULATION_H
#define __SIMULATION_H

#include <vector>
#include <random>
#include "Vec.h"
#include "EventManager.h"

//Start by creating a simple simulation class so we can start testing.
//After sufficient testing, we can move to the entity-component system
//and check how performance is affected. If the performance penalty is
//not too significant, it will beneficial to keep the EC system.

//NOTE: The plan is, for a particle i, to only schedule the next collision
//and eventually, when the cell list is implemented, the next cell transfer.
//check 'Algorithms for Particle-Field Simulations with Collisions'

//NOTE: Consider using smart pointers for passing events around

//LOG: Next time, we want to implement a rayAABBIntersection function for
//calculating when a sphere will cross the simulation box.

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
    bool raySphereIntersection(double radius, const Vec3d& pos, const Vec3d& dir, double& t)const;
    CollisionEvent* getCollisionEvent(size_t pA, size_t pB)const;
    void runCollisionEvent(const CollisionEvent& event);
    Vec3d applyPeriodicBC(const Vec3d& vec)const;
    void updateParticle(size_t pID);

    size_t nSpheres_;
    double boxSize_;
    Time   time_;
    std::vector<size_t> nCollisions_;
    std::vector<double> radii_;
    std::vector<Time>   times_;
    std::vector<Vec3d>  positions_;
    std::vector<Vec3d>  velocities_;

    EventManager eventManager_;

    std::mt19937 mtGen_;
};

#endif
