#include <cstdio>
#include <cmath>
#include "Simulation.h"
#include "shape/variant.h"
#include "shape/box.h"
#include "overlap/ray_casting.h"
#include "overlap/gjk.h"
#include "overlap/overlap.h"
#include "overlap/obb.h"
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

    class BoundingBoxVisitor: public boost::static_visitor<BBShape> {
    public:
        template<typename T>
        BBShape operator()(const T& shape)const{
            double max_x = shape.support(clam::Vec3d(1.0, 0.0, 0.0))[0];
            double min_x = shape.support(clam::Vec3d(-1.0, 0.0, 0.0))[0];
            double max_y = shape.support(clam::Vec3d(0.0, 1.0, 0.0))[1];
            double min_y = shape.support(clam::Vec3d(0.0, -1.0, 0.0))[1];
            double max_z = shape.support(clam::Vec3d(0.0, 0.0, 1.0))[2];
            double min_z = shape.support(clam::Vec3d(0.0, 0.0, -1.0))[2];
            clam::Vec3d max_r = clam::Vec3d(max_x, max_y, max_z);
            clam::Vec3d min_r = clam::Vec3d(min_x, min_y, min_z);

            return BBShape{
                0.5 * (max_r + min_r),
                0.5 * (max_r - min_r)
            };
        }
    };

    BBShape calc_bounding_box(const shape::Variant& shape){
        return boost::apply_visitor(BoundingBoxVisitor(), shape);
    }
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

class Simulation::ShapeBoxCrossEventVisitor: public boost::static_visitor<ParticleEvent> {
public:
    ShapeBoxCrossEventVisitor(const Simulation& sim, int pid):
        sim_(sim), pid_(pid)
    {}

    template<typename T>
    ParticleEvent operator()(const T& shape)const{

        double time(0.0);

        BoundingBox bbox = sim_.boxes_[pid_];
        auto bbox_inv_rot = bbox.rot_.inv();

        Particle part = sim_.particles_[pid_];
        stream_rotation(part, sim_.time_);
        part.pos     = bbox_inv_rot.rotate(sim_.pbc_.minImage(part.pos - bbox.pos_));
        part.rot     = bbox_inv_rot * part.rot;
        part.vel     = bbox_inv_rot.rotate(part.vel);
        part.ang_vel = bbox_inv_rot.rotate(part.ang_vel);
        part.time    = 0.0;

        clam::Vec3d half_size_ = sim_.box_shapes_[part.shape_id].half_size_;
        double out_radius = part.size * shape_outradius(shape);

        //TODO: Look into high iteration number!!!
        while(true){
            double max_advance = std::numeric_limits<double>::max();
            double distance = 0.0;
            int max_id = 0;
            auto inv_rot = part.rot.inv();
            for(int i = 0; i < 3; ++i){
                double max_vel = part.vel[i] + part.ang_vel.length() * out_radius;
                if(max_vel > 0.0){
                    auto dir = clam::Vec3d(0.0);
                    dir[i] = 1.0;
                    auto dir_p = inv_rot.rotate(dir);
                    double dist = half_size_[i] - part.rot.rotate(shape.support(dir_p))[i] - part.pos[i];
                    double advance = dist / max_vel;
                    if(advance < max_advance){
                        distance = dist;
                        max_advance = advance;
                        max_id = 2 * i;
                    }
                }
                max_vel = -part.vel[i] + part.ang_vel.length() * out_radius;
                if(max_vel > 0.0){
                    auto dir = clam::Vec3d(0.0);
                    dir[i] = -1.0;
                    auto dir_p = inv_rot.rotate(dir);
                    double dist = half_size_[i] + part.rot.rotate(shape.support(dir_p))[i] + part.pos[i];
                    double advance = dist / max_vel;
                    if(advance < max_advance){
                        distance = dist;
                        max_advance = advance;
                        max_id = 2 * i + 1;
                    }
                }
            }

            if(distance < sim_.closest_distance_tol_){
                while(distance < sim_.closest_distance_tol_){
                    time *= 0.999;
                    stream_position(part, time);
                    stream_rotation(part, time);
                    part.time = time;

                    auto dir = clam::Vec3d(0.0);
                    int idx = max_id / 2;
                    dir[idx] = 1.0 - 2.0 * (max_id % 2);
                    auto dir_p = part.rot.inv().rotate(dir);
                    distance = half_size_[idx] - dir[idx] * (part.rot.rotate(shape.support(dir_p))[idx] + part.pos[idx]);
                }

                return ParticleEvent::CellCross(sim_.time_ + time, pid_, 0);
            }

            time += max_advance;

            stream_position(part, time);
            stream_rotation(part, time);
            part.time = time;
        }

        return ParticleEvent::None();
    }

    //TODO: Implement
    ParticleEvent operator()(const shape::Sphere& a, const shape::Sphere& b)const{
        return ParticleEvent::None();
    }

private:
    const Simulation& sim_;
    int pid_;
};

ParticleEvent Simulation::get_collision_event(int pA, int pB)const{
    return boost::apply_visitor(
        ShapeCollisionEventVisitor(*this, pA, pB),
        *shapes_[particles_[pA].shape_id], *shapes_[particles_[pB].shape_id]
    );
}

ParticleEvent Simulation::get_cell_cross_event(int pid)const{
    return boost::apply_visitor(ShapeBoxCrossEventVisitor(*this, pid), *shapes_[particles_[pid].shape_id]);
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
    //NOTE: In principle n != pA
    for(int n: nnl_[pA]){
        if(n != pA && n != pB){
            //update_particle(particles_[n], time_, pbc_);
            auto event = get_collision_event(pA, n);
            if(event.get_type() != PE_NONE) event_mgr_.push(pA, event);
        }
    }

    for(int n: nnl_[pB]){
        if(n != pA && n != pB){
            //update_particle(particles_[n], time_, pbc_);
            auto event = get_collision_event(pB, n);
            if(event.get_type() != PE_NONE) event_mgr_.push(pB, event);
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
    int pid = event.pid_;

    for(int nb: nnl_[pid]){
        if(nb == pid) continue;
        //Perhaps remove_if
        nnl_[nb].erase(std::find(nnl_[nb].begin(), nnl_[nb].end(), pid));
    }
    nnl_[pid].clear();

    update_particle(particles_[pid], time_, pbc_);
    boxes_[pid].pos_ = particles_[pid].pos;
    boxes_[pid].rot_ = particles_[pid].rot;

    cll_.update(pid, boxes_[pid].pos_);

    auto shape_id_a = particles_[pid].shape_id;

    for(int cid: cll_.cell_nbs(cll_.cell_index(pid))){
        for(int n: cll_.cell_content(cid)){
            if(n == pid) continue;
            auto shape_id_b = particles_[n].shape_id;
            auto bba = boxes_[pid];
            auto bbb = boxes_[n];
            bbb.pos_ = pbc_.minImage(bbb.pos_ - bba.pos_);
            bba.pos_ = 0.0;
            if(overlap::obb_overlap(bba, box_shapes_[shape_id_a], bbb, box_shapes_[shape_id_b], obb_margin_)){
                nnl_[pid].push_back(n);
                nnl_[n].push_back(pid);
                auto event = get_collision_event(pid, n);
                if(event.get_type() != PE_NONE) event_mgr_.push(pid, event);
            }
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
    obb_margin_(0.1),
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

    //Initialize obb shapes
    box_shapes_ = new BBShape[shapes_.size()];
    for(size_t i = 0; i < shapes_.size(); ++i) box_shapes_[i] = calc_bounding_box(*shapes_[i]);

    boxes_ = new BoundingBox[n_part];
    double max_radius = 0.0;
    //Initialize paricle velocities
    std::uniform_real_distribution<double> dist(-1.0, 1.0);
    clam::Vec3d sys_vel;
    for(size_t i = 0; i < n_part; ++i){
        auto& particle = particles_[i];

        double x1, x2, r;
        do{
            x1 = dist(mtGen_);
            x2 = dist(mtGen_);
            r = x1 * x1 + x2 * x2;
        }while(r >= 1.0);
        double s = 2.0 * sqrt(1.0 - r);

        clam::Vec3d vec(x1 * s, x2 * s, 1.0 - 2.0 * r);
        sys_vel += vec;
        particle.vel = vec;

        particle.pos = pbc_.apply(particle.pos + 0.1); //Make sure pbc are correct at start.
        boxes_[i].pos_ = particle.pos;
        boxes_[i].rot_ = particle.rot;

        //NOTE: Can cache radii for each shape
        double outradius = particle.size * box_shapes_[particle.shape_id].out_radius();
        max_radius = std::max(max_radius, outradius);
    }
    sys_vel = sys_vel * (1.0 / n_part);

    //Initialize cell list
    cll_.init(n_part, pbc_.getSize(), 2.0 * max_radius + 0.01);
    for(size_t i = 0; i < n_part; ++i){
        particles_[i].vel -= sys_vel;
        base_kinetic_energy_ += 0.5 * particles_[i].vel.length2();
        cll_.add(i, particles_[i].pos);
    }

    //Initialize nnl
    nnl_ = new std::vector<size_t>[n_part];

    foreach_pair(cll_, [this](int i, int j) -> bool {
        auto shape_id_i = particles_[i].shape_id;
        auto shape_id_j = particles_[j].shape_id;
        auto bbi = boxes_[i];
        auto bbj = boxes_[j];
        bbj.pos_ = pbc_.minImage(boxes_[j].pos_ - boxes_[i].pos_);
        bbi.pos_ = 0.0;
        if(overlap::obb_overlap(bbi, box_shapes_[shape_id_i], bbj, box_shapes_[shape_id_j], obb_margin_)){
            nnl_[i].push_back(j);
            nnl_[j].push_back(i);
            auto event = get_collision_event(i, j);
            if(event.get_type() != PE_NONE) event_mgr_.push(i, event);
        }
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
        assert(event.get_type() != PE_NONE);
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
