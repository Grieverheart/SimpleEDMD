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
#include "bounding_volume_variant_fwd.h"
#include "serialization/archive.h"

namespace shape{
    class Box;
}

struct Transform;

class Simulation{
public:
    Simulation(const Configuration& config);
    Simulation(void);
    ~Simulation(void);

    void run(double end_time, PeriodicCallback& output_condition);
    void stop(void);
    void restart(void);

    const std::vector<Particle>& particles(void)const;
    const RectangularPBC& pbc(void)const;
    const Configuration& configuration(void)const;

    void reset_statistics(void);

    double time(void)const;
    double average_stress(void)const;
    double average_kinetic_energy(void)const;
    double average_pressure(void)const;
    int num_collisions(void)const;

    bool check_overlaps(void)const;

    friend void serialize(Archive&, const Simulation&);
    friend void deserialize(Archive&, Simulation*);

private:
    ParticleEvent get_collision_event(int pA, int pB)const;
    ParticleEvent get_neighborhood_cross_event(int pid)const;
    void run_collision_event(const ParticleEvent& event);
    void run_neighborhood_cross_event(const ParticleEvent& event);
    void run_possible_collision_event(const ParticleEvent& event);

private:
    bool is_running_;
    //bool was_stopped_;

    double time_;
    double prev_time_;
    double statistics_start_time_;
    double closest_distance_tol_;
    double obb_margin_;
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

    //NOTE: We could make this into a separate class
    //but not sure how to handle shape_id without
    //duplicating data.
    bounding_volume::Variant** box_shapes_;
    Transform* boxes_;
    std::vector<size_t>* nnl_;

    EventManager event_mgr_;
    CellList     cll_;

    //TODO: Replace with original mt19937
    std::mt19937 mtGen_;

    class ShapeCollisionEventVisitor;
    class ShapeBoxCrossEventVisitor;
};

#endif
