#include <cstdio>
#include <cmath>
#include "Simulation.h"
#include "shape/variant.h"
#include "overlap/ray_casting.h"
#include "overlap/gjk.h"
#include "overlap/overlap.h"
#include "serialization/common.h"
#include "serialization/vector.h"

static void print_particle(const Particle& p){
    printf("pos: glm::vec3(%e, %e, %e)\n", p.pos[0], p.pos[1], p.pos[2]);
    clam::Vec3d axis;
    double angle;
    p.rot.toAxisAngle(angle, axis);
    printf("rot: %e, glm::vec3(%e, %e, %e)\n", angle, axis[0], axis[1], axis[2]);
    printf("vel: glm::vec3(%e, %e, %e)\n", p.vel[0], p.vel[1], p.vel[2]);
    printf("ang_vel: glm::vec3(%e, %e, %e)\n", p.ang_vel[0], p.ang_vel[1], p.ang_vel[2]);
}

namespace{
    template<typename F>
    inline bool foreach_pair(const CellList& cll_, F func){
        for(auto cell: cll_.cells()){
            for(auto cell_nb: cll_.cell_vol_nbs(cell)){
                if(cell == cell_nb){
                    for(auto pa_itr = cll_.cell_content(cell); pa_itr != pa_itr.end(); ++pa_itr){
                        for(auto pb_itr = pa_itr + 1; pb_itr != pb_itr.end(); ++pb_itr){
                            if(func(*pa_itr, *pb_itr)) return true;
                        }
                    }
                }
                else{
                    for(auto pa: cll_.cell_content(cell)){
                        for(auto pb: cll_.cell_content(cell_nb)){
                            if(func(pa, pb)) return true;
                        }
                    }
                }
            }
        }
        return false;
    }

    //TODO: Move these to a utility file so they can be used by others.
    inline void stream_position(Particle& particle, double time){
        particle.pos += particle.vel * (time - particle.time);
    }

    inline void stream_rotation(Particle& particle, double time){
        clam::Vec3d ha = ((time - particle.time) * 0.5) * particle.ang_vel;
        double l = ha.length(); // magnitude
        if(l > 0.0){
            double sl, cl;
            sincos(l, &sl, &cl);
            particle.rot = clam::Quatd(ha * (sl / l), cl) * particle.rot;
        }
    }

    inline void update_particle(Particle& particle, double time, const RectangularPBC& pbc){
        if(particle.time < time){
            stream_position(particle, time);
            stream_rotation(particle, time);
            particle.pos = pbc.apply(particle.pos);
            particle.time = time;
        }
    }

    inline double max_overlap_distance(double size_a, const shape::Variant& sa, double size_b, const shape::Variant& sb){
        return size_a * shape_outradius(sa) + size_b * shape_outradius(sb);
    }

    class ShapeCollisionResolutionVisitor: public boost::static_visitor<bool> {
    public:
        ShapeCollisionResolutionVisitor(Particle& pa, Particle& pb, const RectangularPBC& pbc, double& momentum_transfer, double& kinetic_delta):
            pa_(pa), pb_(pb), pbc_(pbc),
            momentum_transfer_(momentum_transfer),
            kinetic_delta_(kinetic_delta)
        {}

        template<typename T, typename U>
        bool operator()(const T& a, const U& b)const{
            clam::Vec3d rel_pos = pbc_.minImage(pb_.pos - pa_.pos);
            Particle temp_pa = pa_;
            Particle temp_pb = pb_;
            temp_pb.pos = rel_pos;
            temp_pa.pos = 0.0;

            clam::Vec3d point_on_a, point_on_b;
            clam::Vec3d dist_vec = overlap::gjk_closest_points(temp_pa, a, temp_pb, b, point_on_a, point_on_b);
            if(dist_vec.length2() == 0.0) return false;

            //@note: We're using the distance vector itself for better accuracy.
            clam::Vec3d normal = -dist_vec / dist_vec.length();

            point_on_b -= rel_pos;
            clam::Vec3d rel_vel = pb_.vel - pa_.vel + clam::cross(pb_.ang_vel, point_on_b) - clam::cross(pa_.ang_vel, point_on_a);
            clam::Vec3d tangent_a = clam::cross(point_on_a, normal);
            clam::Vec3d tangent_b = clam::cross(point_on_b, normal);

            double momentum_delta =
                2.0 * clam::dot(normal, rel_vel) /
                (2.0 + tangent_a.length2() + tangent_b.length2());

            momentum_transfer_ = momentum_delta * clam::dot(normal, rel_pos);

            clam::Vec3d vec_momentum_delta = momentum_delta * normal;
            kinetic_delta_ = clam::dot(vec_momentum_delta, vec_momentum_delta + pa_.vel - pb_.vel);

            pa_.vel += vec_momentum_delta;
            pb_.vel -= vec_momentum_delta;

            pa_.ang_vel += momentum_delta * tangent_a;
            pb_.ang_vel -= momentum_delta * tangent_b;

            return true;
        }

        //Assume always resolved.
        //TODO: Handle momentum_transfer
        bool operator()(const shape::Sphere& a, const shape::Sphere& b)const{
            clam::Vec3d rel_pos = pbc_.minImage(pa_.pos - pb_.pos);
            double dist2 = rel_pos.length2();

            clam::Vec3d rel_vel = pa_.vel - pb_.vel;
            clam::Vec3d delta_vel = rel_pos * clam::dot(rel_pos, rel_vel) / dist2;

            pa_.vel -= delta_vel;
            pb_.vel += delta_vel;

            return true;
        }

    private:
        Particle& pa_;
        Particle& pb_;
        const RectangularPBC& pbc_;
        double& momentum_transfer_;
        double& kinetic_delta_;
    };
}

class Simulation::ShapeCollisionEventVisitor: public boost::static_visitor<ParticleEvent> {
public:
    ShapeCollisionEventVisitor(const Simulation& sim, int pa_idx, int pb_idx):
        sim_(sim), pa_idx_(pa_idx), pb_idx_(pb_idx)
    {}

    template<typename T, typename U>
    ParticleEvent operator()(const T& a, const U& b)const{
        Particle partA = sim_.particles_[pa_idx_];
        Particle partB = sim_.particles_[pb_idx_];
        clam::Vec3d dist = sim_.pbc_.minImage(partB.pos + partB.vel * (sim_.time_ - partB.time) - partA.pos);
        double max_dist = max_overlap_distance(partA.size, *sim_.shapes_[partA.shape_id], partB.size, *sim_.shapes_[partB.shape_id]);

        if(dist.length() > max_dist + 2.0 * sim_.closest_distance_tol_){
            double time(0.0);
            clam::Vec3d rel_vel = partA.vel - partB.vel;
            if(overlap::sphere_raycast(max_dist + 1.5 * sim_.closest_distance_tol_, dist, rel_vel, time)){
                return ParticleEvent::PossibleCollision(sim_.time_ + time, pa_idx_, pb_idx_, sim_.n_collisions_[pb_idx_]);
            }
        }
        else{
            double time(0.0);

            stream_rotation(partB, sim_.time_);
            partB.pos  = dist;
            partB.vel  = partB.vel - partA.vel;
            partB.time = 0.0;

            partA.pos  = 0.0;
            partA.vel  = 0.0;
            partA.time = 0.0;

            //TODO: Look into high iteration number!!!
            while(true){
                clam::Vec3d shortest_dist = overlap::gjk_distance(partA, a, partB, b, sim_.closest_distance_tol_);

                double distance = shortest_dist.length();
                if(distance == 0.0){
                    shortest_dist = overlap::gjk_distance(partA, a, partB, b);
                    distance = shortest_dist.length();
                }
                clam::Vec3d shortest_dist_n = shortest_dist / distance;

                /* Conservative advancement by Mirtich 1996 PhD Thesis */
                //double out_radius_A = partA.size * shape_outradius(a);
                //double out_radius_B = partB.size * shape_outradius(b);
                //double max_vel = clam::dot(shortest_dist_n, partB.vel) +
                //                 partA.ang_vel.length() * out_radius_A + partB.ang_vel.length() * out_radius_B;


                /* Conservative advancement, tighter bound using support functions. */
                //NOTE: Most probably is only correct for point symmetric particle shapes.
                double max_vel = clam::dot(shortest_dist_n, partB.vel);
                double abs_omega_A = partA.ang_vel.length();
                if(abs_omega_A > 0.0){
                    auto inv_rot_A = partA.rot.inv();
                    clam::Vec3d c1  = inv_rot_A.rotate(clam::cross(partA.ang_vel, shortest_dist_n));
                    clam::Vec3d c1p = clam::cross(inv_rot_A.rotate(partA.ang_vel), c1) / abs_omega_A;
                    max_vel += clam::dot(c1, partA.size * a.support(c1)) + clam::dot(c1p, partA.size * a.support(c1p));
                }
                double abs_omega_B = partB.ang_vel.length();
                if(abs_omega_B > 0.0){
                    auto inv_rot_B = partB.rot.inv();
                    clam::Vec3d c2  = inv_rot_B.rotate(clam::cross(partB.ang_vel, shortest_dist_n));
                    clam::Vec3d c2p = clam::cross(inv_rot_B.rotate(partB.ang_vel), c2) / abs_omega_B;
                    max_vel += clam::dot(c2, partB.size * b.support(c2)) + clam::dot(c2p, partB.size * b.support(c2p));
                }

                if(max_vel < 0.0) break;

                //TODO: We can change the bound to 2 * tol and additionally check
                //if we are approaching (distance < prev_distance).
                if(distance < sim_.closest_distance_tol_){
                    int iter = 0;
                    while(shortest_dist.length() < sim_.closest_distance_tol_){
                        time *= 0.999;
                        stream_position(partB, time);
                        stream_rotation(partB, time);
                        stream_rotation(partA, time);
                        partB.time = time;
                        partA.time = time;
                        shortest_dist = overlap::gjk_distance(partA, a, partB, b, sim_.closest_distance_tol_);
                        //NOTE: This should alsmost never happen.
                        if(iter++ > 1000) return ParticleEvent::None();
                    }

                    return ParticleEvent::Collision(sim_.time_ + time, pa_idx_, pb_idx_, sim_.n_collisions_[pb_idx_]);
                }

                double max_advance = distance / max_vel;
                time += max_advance;

                //NOTE: The separation time of the bounding spheres is a
                //convenient but too conservative bound for the upper bound
                //to the collision time.
                if(time < 0.0 || time > sim_.max_collision_time_) break;

                stream_position(partB, time);
                stream_rotation(partB, time);
                stream_rotation(partA, time);
                partB.time = time;
                partA.time = time;
            }
        }

        return ParticleEvent::None();
    }

    ParticleEvent operator()(const shape::Sphere& a, const shape::Sphere& b)const{
        const Particle& partA = sim_.particles_[pa_idx_];
        const Particle& partB = sim_.particles_[pb_idx_];

        clam::Vec3d dist = sim_.pbc_.minImage(partB.pos + partB.vel * (sim_.time_ - partB.time) - partA.pos);
        clam::Vec3d rel_vel = partA.vel - partB.vel;

        double time(0.0);
        clam::Vec3d normal;
        if(overlap::sphere_raycast(partA.size * a.radius() + partB.size * b.radius() + sim_.closest_distance_tol_, dist, rel_vel, time, &normal)){
            return ParticleEvent::Collision(sim_.time_ + time, pa_idx_, pb_idx_, sim_.n_collisions_[pb_idx_]);
        }
        else return ParticleEvent::None();
    }

private:
    const Simulation& sim_;
    int pa_idx_;
    int pb_idx_;
};

ParticleEvent Simulation::get_collision_event(int pA, int pB)const{
    return boost::apply_visitor(
        ShapeCollisionEventVisitor(*this, pA, pB),
        *shapes_[particles_[pA].shape_id], *shapes_[particles_[pB].shape_id]
    );
}

ParticleEvent Simulation::get_cell_cross_event(int pid)const{
    int cidx = cll_.cell_index(pid);
    double time(0.0);
    clam::Vec3d rpos = pbc_.minImage(particles_[pid].pos - cll_.cell_origin(cidx));
    int cellOffset = overlap::cell_raycast(cll_.cell_size(), rpos, particles_[pid].vel, time);
    return ParticleEvent::CellCross(time_ + time, pid, cellOffset);
}

//NOTE: For simplicity, for now we assume equal mass spheres
void Simulation::run_collision_event(const ParticleEvent& event){
    int pA = event.pid_;
    int pB = event.get_id();

    if(n_collisions_[pB] != event.optional_){
        event_mgr_.update(pA);
        return;
    }

    {
        double time_diffA = time_ - particles_[pA].time;
        double time_diffB = time_ - particles_[pB].time;
        if(time_diffA > max_inflight_time_) max_inflight_time_ = time_diffA;
        if(time_diffB > max_inflight_time_) max_inflight_time_ = time_diffB;
    }

    update_particle(particles_[pA], time_, pbc_);
    update_particle(particles_[pB], time_, pbc_);

    //Resolve collision.
    double momentum_transfer = 0.0;
    double kinetic_delta = 0.0;
    bool resolved = boost::apply_visitor(
        ShapeCollisionResolutionVisitor(particles_[pA], particles_[pB], pbc_, momentum_transfer, kinetic_delta),
        *shapes_[particles_[pA].shape_id], *shapes_[particles_[pB].shape_id]
    );

    if(!resolved){
        double temp_time = (time_ + prev_time_) * 0.5;
        stream_position(particles_[pA], temp_time);
        stream_position(particles_[pB], temp_time);
        stream_rotation(particles_[pA], temp_time);
        stream_rotation(particles_[pB], temp_time);
        particles_[pA].time = temp_time;
        particles_[pB].time = temp_time;
        time_ = temp_time;
        auto event = get_collision_event(pA, pB);
        if(event.get_type() != PE_NONE){
            event_mgr_.push(pA, event);
            event_mgr_.update(pA);
        }
        return;
    }

    av_momentum_transfer_ += momentum_transfer;
    kinetic_delta_ += kinetic_delta;

    ++n_collisions_[pA];
    ++n_collisions_[pB];

    event_mgr_.clear(pA);
    event_mgr_.clear(pB);

    {
        auto event = get_collision_event(pA, pB);
        if(event.get_type() != PE_NONE){
            if(event.time_ > time_) event_mgr_.push(pA, event);
            else assert(false);
        }
    }

    //Recalculate collision events
    for(int cid: cll_.particle_cell_nbs(pA)){
        for(int n: cll_.cell_content(cid)){
            if(n != pA && n != pB){
                //update_particle(particles_[n], time_, pbc_);
                auto event = get_collision_event(pA, n);
                if(event.get_type() != PE_NONE) event_mgr_.push(pA, event);
            }
        }
    }

    for(int cid: cll_.particle_cell_nbs(pB)){
        for(int n: cll_.cell_content(cid)){
            if(n != pA && n != pB){
                //update_particle(particles_[n], time_, pbc_);
                auto event = get_collision_event(pB, n);
                if(event.get_type() != PE_NONE) event_mgr_.push(pB, event);
            }
        }
    }

    event_mgr_.push(pA, get_cell_cross_event(pA));
    event_mgr_.push(pB, get_cell_cross_event(pB));

    event_mgr_.update(pA);
    event_mgr_.update(pB);
}

void Simulation::run_possible_collision_event(const ParticleEvent& event){
    int pA = event.pid_;
    int pB = event.get_id();

    if(n_collisions_[pB] != event.optional_){
        event_mgr_.update(pA);
        return;
    }

    update_particle(particles_[pA], time_, pbc_);

    auto new_event = get_collision_event(pA, pB);
    if(new_event.get_type() != PE_NONE){
        event_mgr_.push(pA, new_event);
    }

    event_mgr_.update(pA);
}

void Simulation::run_cell_cross_event(const ParticleEvent& event){
    int pid     = event.pid_;
    int coffset = event.get_id();
    cll_.move(pid, coffset);
    update_particle(particles_[pid], time_, pbc_);
    for(int cid: cll_.cell_dir_nbs(cll_.cell_index(pid), coffset)){
        for(int n: cll_.cell_content(cid)){
            //update_particle(particles_[n], time_, pbc_);
            auto event = get_collision_event(pid, n);
            if(event.get_type() != PE_NONE) event_mgr_.push(pid, event);
        }
    }
    event_mgr_.push(pid, get_cell_cross_event(pid));
    event_mgr_.update(pid);
}

const std::vector<Particle>& Simulation::particles(void)const{
    return particles_;
}

const RectangularPBC& Simulation::pbc(void)const{
    return pbc_;
}

const Configuration& Simulation::configuration(void)const{
    return config_;
}
Simulation::Simulation(void):
    pbc_(config_.pbc_), particles_(config_.particles_), shapes_(config_.shapes_)
{}

//TODO: Perhaps add exceptions to constructor
Simulation::Simulation(const Configuration& config):
    time_(0.0),
    statistics_start_time_(0.0),
    closest_distance_tol_(1.0e-8), //@note: increase tolerance to increase performance.
    max_collision_time_(5.0),
    av_momentum_transfer_(0.0),
    av_kinetic_delta_(0.0), base_kinetic_energy_(0.0), kinetic_delta_(0.0),
    max_inflight_time_(0.0),
    n_collision_events_(0),
    config_(config),
    pbc_(config_.pbc_), particles_(config_.particles_), shapes_(config_.shapes_)
{
    mtGen_.seed(0);//time(NULL));

    auto n_part = particles_.size();

    if(n_part) event_mgr_.resize(n_part);
    else return;

    //Initialize number of collisions to zero
    n_collisions_.resize(n_part, 0);

    //Initialize cell list
    double max_radius = 0.0;
    for(auto particle: particles_){
        double outradius = particle.size * boost::apply_visitor(ShapeOutRadiusVisitor(), *shapes_[particle.shape_id]);
        max_radius = std::max(max_radius, outradius);
    }
    cll_.init(n_part, pbc_.getSize(), 2.0 * max_radius + 0.01);

    //Initialize paricle velocities
    std::uniform_real_distribution<double> dist(-1.0, 1.0);
    clam::Vec3d sys_vel;
    for(size_t i = 0; i < n_part; ++i){
        double x1, x2, r;
        do{
            x1 = dist(mtGen_);
            x2 = dist(mtGen_);
            r = x1 * x1 + x2 * x2;
        }while(r >= 1.0);
        double s = 2.0 * sqrt(1.0 - r);
        clam::Vec3d vec(x1 * s, x2 * s, 1.0 - 2.0 * r);

        //IMPORTANT: If the particle lies exactly on the cell boundary, we get
        //errors in the simulation. The time for the next cell crossing is
        //correctly found to be 0.0, but the event calendar cannot handle
        //simultaneous events.
        sys_vel += vec;
        particles_[i].vel = vec;
        particles_[i].pos = pbc_.apply(particles_[i].pos + 0.1); //Make sure pbc are correct at start.
        cll_.add(i, particles_[i].pos);
    }
    sys_vel = sys_vel * (1.0 / n_part);

    for(size_t i = 0; i < n_part; ++i){
        particles_[i].vel -= sys_vel;
        base_kinetic_energy_ += 0.5 * particles_[i].vel.length2();
    }

    foreach_pair(cll_, [this](int i, int j) -> bool {
        auto event = get_collision_event(i, j);
        if(event.get_type() != PE_NONE) event_mgr_.push(i, event);
//Check for overlaps at start of simulation.
#ifndef NDEBUG
        Particle pa = particles_[i];
        Particle pb = particles_[j];
        pb.pos = pbc_.minImage(pb.pos - pa.pos);
        pa.pos = 0.0;
        assert(overlap::shape_overlap(pa, *shapes_[pa.shape_id], pb, *shapes_[pb.shape_id]) == false);
#endif
        return false;
    });

    for(size_t i = 0; i < n_part; ++i){
        auto event = get_cell_cross_event(i);
        event_mgr_.push(i, event);
    }

    event_mgr_.init();
}

void Simulation::run(double end_time, PeriodicCallback& output_condition){

    bool running = true;
    while(running){
        ParticleEvent next_event = event_mgr_.getNextEvent();
        prev_time_ = time_;
        time_ = next_event.time_;

        av_kinetic_delta_ += (time_ - prev_time_) * kinetic_delta_;

        switch(next_event.get_type()){
        case PE_COLLISION:
            ++n_collision_events_;
            run_collision_event(next_event);
            break;
        case PE_POSSIBLE_COLLISION:
            run_possible_collision_event(next_event);
            break;
        case PE_CELLCROSS:
            run_cell_cross_event(next_event);
            break;
        default:
            running = false;
            break;
        }

        output_condition(time_);

        if(time_ >= end_time) running = false;

        //Optimization: After we have collected enough events, we set the
        //maximum time for which to search for collisions to 1.5 times
        //the maximum in-flight time of each particle. Of course, it
        //might still be possible to miss a collision, but unlikely.
        if(n_collision_events_ == 40000) max_collision_time_ = 1.5 * max_inflight_time_;
    }
}

void Simulation::reset_statistics(void){
    statistics_start_time_ = time_;
    av_kinetic_delta_      = 0.0;
    av_momentum_transfer_  = 0.0;
    n_collision_events_    = 0;
}

int Simulation::num_collisions(void)const{
    return n_collision_events_;
}

double Simulation::time(void)const{
    return time_;
}

double Simulation::average_stress(void)const{
    return av_momentum_transfer_ / (time_ - statistics_start_time_);
}

double Simulation::average_kinetic_energy(void)const{
    return base_kinetic_energy_ + av_kinetic_delta_ / (time_ - statistics_start_time_);
}

double Simulation::average_pressure(void)const{
    clam::Vec3d box_size = pbc_.getSize();
    double volume = box_size[0] * box_size[1] * box_size[2];
    double kT = 2.0 * average_kinetic_energy() / (3.0 * particles_.size());
    double pressure = (particles_.size() - average_stress() / (3.0 * kT)) / volume;
    return pressure;
}

void serialize(Archive& ar, const Simulation& sim){
    serialize(ar, sim.time_);
    serialize(ar, sim.prev_time_);
    serialize(ar, sim.statistics_start_time_);
    serialize(ar, sim.closest_distance_tol_);
    serialize(ar, sim.max_collision_time_);
    serialize(ar, sim.av_momentum_transfer_);
    serialize(ar, sim.av_kinetic_delta_);
    serialize(ar, sim.base_kinetic_energy_);
    serialize(ar, sim.kinetic_delta_);
    serialize(ar, sim.max_inflight_time_);
    serialize(ar, sim.n_collision_events_);
    serialize(ar, sim.n_collisions_);
    serialize(ar, sim.config_);
    serialize(ar, sim.event_mgr_);
    serialize(ar, sim.cll_);
}

void deserialize(Archive& ar, Simulation* sim){
    deserialize(ar, &sim->time_);
    deserialize(ar, &sim->prev_time_);
    deserialize(ar, &sim->statistics_start_time_);
    deserialize(ar, &sim->closest_distance_tol_);
    deserialize(ar, &sim->max_collision_time_);
    deserialize(ar, &sim->av_momentum_transfer_);
    deserialize(ar, &sim->av_kinetic_delta_);
    deserialize(ar, &sim->base_kinetic_energy_);
    deserialize(ar, &sim->kinetic_delta_);
    deserialize(ar, &sim->max_inflight_time_);
    deserialize(ar, &sim->n_collision_events_);
    deserialize(ar, &sim->n_collisions_);
    deserialize(ar, &sim->config_);
    deserialize(ar, &sim->event_mgr_);
    deserialize(ar, &sim->cll_);
}
