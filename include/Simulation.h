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
    Simulation(const Configuration& config);

    void run(double end_time, PeriodicCallback& output_condition);

    const std::vector<Particle>& get_particles(void)const;
    const CubicPBC& get_pbc(void)const;
    clam::Vec3d get_system_velocity(void)const;
    const Configuration& get_configuration(void)const;

private:
    ParticleEvent get_collision_event(int pA, int pB)const;
    ParticleEvent get_cell_cross_event(int pid)const;
    void run_collision_event(const ParticleEvent& event);
    void run_cell_cross_event(const ParticleEvent& event);
    void run_possible_collision_event(const ParticleEvent& event);

private:
    double time_;
    double prev_time_;
    double closest_distance_tol_;
    clam::Vec3d sys_vel_;
    std::vector<uint32_t> n_collisions_;

    Configuration config_;
    CubicPBC& pbc_;
    std::vector<Particle>& particles_;
    std::vector<shape::Variant*>& shapes_;

    EventManager event_mgr_;
    CellList     cll_;

    std::mt19937 mtGen_;

    class ShapeCollisionEventVisitor;
};

#endif
