#ifndef __SIMULATION_H
#define __SIMULATION_H

#include <vector>
#include <random>
#include <ctime>
#include "BoundaryCondition.h"
#include "PeriodicCallback.h"
#include "CellList.h"
#include "particle.h"
#include "EventManager.h"
#include "shape/variant_fwd.h"

//NOTE:
//1. The Boundary Condition class will most probably need to be handed
//   as a template parameter for performance reasons.
//2. Make a Collision Response class for different collision responses.
//   We might have to pass this as a template parameter too.
//3. We might also need to switch to -B/2 -> +B/2 coordinates so that
//   we can easily implement spherical boundary conditions.

class TempEventManager;
class EventManager;

class Simulation{
public:
    Simulation(const CubicPBC& pbc, std::vector<Particle>&& particles, std::vector<shape::Variant*>&& shapes):
        n_part_(particles.size()), n_shapes_(shapes.size()), time_(0.0),
        closest_distance_tol_(1.0e-10), //@note: increase tolerance to increase performance.
        pbc_(pbc), particles_(particles), shapes_(shapes), systemVelocity_(0.0)
    {
        mtGen_.seed(0);//time(NULL));
    }

    void run(double endTime, PeriodicCallback& outputCondition);
    bool init(void);

    int getNumParticles(void)const{
        return n_part_;
    }
    const std::vector<Particle>& getParticles(void)const{
        return particles_;
    }
    const CubicPBC& getPBC(void)const{
        return pbc_;
    }
    clam::Vec3d getSystemVelocity(void)const{
        return systemVelocity_;
    }

private:
    ParticleEvent getCollisionEvent(int pA, int pB)const;
    ParticleEvent getCellCrossEvent(int pid)const;
    void runCollisionEvent(const ParticleEvent& event);
    void runCellCrossEvent(const ParticleEvent& event);
    void runPossibleCollisionEvent(const ParticleEvent& event);
    void updateParticle(int pid);

    int    n_part_;
    int    n_shapes_;
    double time_;
    double prev_time_;
    double closest_distance_tol_;

    CubicPBC pbc_;

    std::vector<uint32_t> nCollisions_;
    std::vector<Particle> particles_;
    std::vector<shape::Variant*> shapes_;
    clam::Vec3d systemVelocity_;

    EventManager eventManager_;
    CellList     cll_;

    std::mt19937 mtGen_;

    class ShapeCollisionEventVisitor;
};

#endif
