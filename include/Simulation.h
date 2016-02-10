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
#include "serialization/archive.h"

class Simulation{
public:
    Simulation(const Configuration& config);
    Simulation(void);

    void run(double end_time, PeriodicCallback& output_condition);

    const std::vector<Particle>& particles(void)const;
    const RectangularPBC& pbc(void)const;
    const Configuration& configuration(void)const;

    void reset_statistics(void);

    double time(void)const;
    double average_stress(void)const;
    double average_kinetic_energy(void)const;
    double average_pressure(void)const;
    int num_collisions(void)const;

    friend void serialize(Archive&, const Simulation&);
    friend void deserialize(Archive&, Simulation*);

private:
    ParticleEvent get_collision_event(int pA, int pB)const;
    ParticleEvent get_cell_cross_event(int pid)const;
    void run_collision_event(const ParticleEvent& event);
    void run_cell_cross_event(const ParticleEvent& event);
    void run_possible_collision_event(const ParticleEvent& event);

private:
    double time_;
    double prev_time_;
    double statistics_start_time_;
    double closest_distance_tol_;
    double max_collision_time_;

    double av_momentum_transfer_;

    //Split kinetic energy into base + delta for better accuracy.
    double av_kinetic_delta_;
    double base_kinetic_energy_;
    double kinetic_delta_;

    double max_inflight_time_;

    int n_collision_events_;

    std::vector<uint32_t> n_collisions_;

    Configuration config_;
    RectangularPBC& pbc_;
    std::vector<Particle>& particles_;
    std::vector<shape::Variant*>& shapes_;

    EventManager event_mgr_;
    CellList     cll_;

    //NOTE: This is not being serialized as it only used at the start of
    //the simulation. In any case, it doesn't hurt if we reseed it.
    std::mt19937 mtGen_;

    class ShapeCollisionEventVisitor;
};

#endif
