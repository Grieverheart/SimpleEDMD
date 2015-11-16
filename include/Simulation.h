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
#include "configuration.h"
#include "shape/variant_fwd.h"

class Simulation{
public:
    Simulation(const Configuration& config):
        time_(0.0),
        closest_distance_tol_(1.0e-10), //@note: increase tolerance to increase performance.
        systemVelocity_(0.0),
        config_(config),
        pbc_(config_.pbc_), particles_(config_.particles_), shapes_(config_.shapes_)
    {
        mtGen_.seed(0);//time(NULL));
    }

    void run(double endTime, PeriodicCallback& outputCondition);
    bool init(void); //TODO: Perhaps move to constructor.

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

private:
    double time_;
    double prev_time_;
    double closest_distance_tol_;
    clam::Vec3d systemVelocity_;
    std::vector<uint32_t> nCollisions_;

    Configuration config_;
    CubicPBC& pbc_;
    std::vector<Particle>& particles_;
    std::vector<shape::Variant*>& shapes_;

    EventManager eventManager_;
    CellList     cll_;

    std::mt19937 mtGen_;

    class ShapeCollisionEventVisitor;
};

#endif
