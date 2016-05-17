#include <cstdio>
#include <cmath>
#include "Simulation.h"
#include "shape/variant.h"
#include "bounding_volume_variant.h"
#include "shape/box.h"
#include "overlap/ray_casting.h"
#include "overlap/gjk.h"
#include "overlap/overlap.h"
#include "overlap/bv_overlap.h"
#include "serialization/common.h"
#include "serialization/vector.h"
#include "serialization/bounding_volume_variant.h"
#include "transform.h"

#define _UNUSED(x) ((void)(x))

static void print_particle(const Particle& p){
    printf("pos: glm::vec3(%e, %e, %e)\n", p.xform.pos_[0], p.xform.pos_[1], p.xform.pos_[2]);
    clam::Vec3d axis;
    double angle;
    p.xform.rot_.toAxisAngle(angle, axis);
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
        particle.xform.pos_ += particle.vel * (time - particle.time);
    }

    inline void stream_rotation(Particle& particle, double time){
        clam::Vec3d ha = ((time - particle.time) * 0.5) * particle.ang_vel;
        double l = ha.length(); // magnitude
        if(l > 0.0){
            double sl = sin(l);
            double cl = cos(l);
            particle.xform.rot_ = clam::Quatd(ha * (sl / l), cl) * particle.xform.rot_;
        }
    }

    inline void update_particle(Particle& particle, double time, const RectangularPBC& pbc){
        if(particle.time < time){
            stream_position(particle, time);
            stream_rotation(particle, time);
            particle.xform.pos_ = pbc.apply(particle.xform.pos_);
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

        //TODO: For spheres we don't need to update angular velocity.
        template<typename T, typename U>
        bool operator()(const T& a, const U& b)const{
            clam::Vec3d rel_pos = pbc_.minImage(pb_.xform.pos_ - pa_.xform.pos_);
            Transform temp_pa = pa_.xform;
            Transform temp_pb = pb_.xform;
            temp_pb.pos_ = rel_pos;
            temp_pa.pos_ = 0.0;

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
            clam::Vec3d rel_pos = pbc_.minImage(pa_.xform.pos_ - pb_.xform.pos_);
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

    class BoundingBoxVisitor: public boost::static_visitor<bounding_volume::Variant> {
    public:
        template<typename T>
        bounding_volume::Variant operator()(const T& shape)const{
            double max_x = shape.support(clam::Vec3d(1.0, 0.0, 0.0))[0];
            double min_x = shape.support(clam::Vec3d(-1.0, 0.0, 0.0))[0];
            double max_y = shape.support(clam::Vec3d(0.0, 1.0, 0.0))[1];
            double min_y = shape.support(clam::Vec3d(0.0, -1.0, 0.0))[1];
            double max_z = shape.support(clam::Vec3d(0.0, 0.0, 1.0))[2];
            double min_z = shape.support(clam::Vec3d(0.0, 0.0, -1.0))[2];
            clam::Vec3d max_r = clam::Vec3d(max_x, max_y, max_z);
            clam::Vec3d min_r = clam::Vec3d(min_x, min_y, min_z);

            return bounding_volume::Variant(shape::Box(clam::Vec3d(
                2.0 * std::max(fabs(min_r[0]), fabs(max_r[0])),
                2.0 * std::max(fabs(min_r[1]), fabs(max_r[1])),
                2.0 * std::max(fabs(min_r[2]), fabs(max_r[2]))
            )));
        }

        bounding_volume::Variant operator()(const shape::Sphere& sphere)const{
            return bounding_volume::Variant(shape::Sphere());
        }
    };

    bounding_volume::Variant calc_bounding_volume(const shape::Variant& shape){
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
        clam::Vec3d dist = sim_.pbc_.minImage(partB.xform.pos_ + partB.vel * (sim_.time_ - partB.time) - partA.xform.pos_);
        double max_dist = max_overlap_distance(partA.xform.size_, *sim_.shapes_[partA.shape_id], partB.xform.size_, *sim_.shapes_[partB.shape_id]);

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
            partB.xform.pos_ = dist;
            partB.vel        = partB.vel - partA.vel;
            partB.time       = 0.0;

            partA.xform.pos_ = 0.0;
            partA.vel        = 0.0;
            partA.time       = 0.0;

            double prev_distance = 0.0;
            int iters = 0;
            while(true){
                ++iters;
                clam::Vec3d shortest_dist = overlap::gjk_distance(partA.xform, a, partB.xform, b);
                double distance = shortest_dist.length();

                if(distance < prev_distance && distance < 2.0 * sim_.closest_distance_tol_){
                    int iter = 0;
                    while(distance < 2.0 * sim_.closest_distance_tol_){
                        time *= 0.999;
                        stream_position(partB, time);
                        stream_rotation(partB, time);
                        stream_rotation(partA, time);
                        partB.time = time;
                        partA.time = time;
                        distance = overlap::gjk_distance(partA.xform, a, partB.xform, b).length();
                        //NOTE: This should alsmost never happen.
                        if(iter++ > 10000){
                            printf("%d, %e, %e, %e\n", iters, prev_distance, distance, time);
                            printf("%d, %d\n", pa_idx_, pb_idx_);
                            return ParticleEvent::None();
                        }
                    }

                    return ParticleEvent::Collision(sim_.time_ + time, pa_idx_, pb_idx_, sim_.n_collisions_[pb_idx_]);
                }

                clam::Vec3d shortest_dist_n = shortest_dist / distance;

                /* Conservative advancement by Mirtich 1996 PhD Thesis */
                //double out_radius_A = partA.xform.size_ * shape_outradius(a);
                //double out_radius_B = partB.xform.size_ * shape_outradius(b);
                //double max_vel = clam::dot(shortest_dist_n, partB.vel) +
                //                 partA.ang_vel.length() * out_radius_A + partB.ang_vel.length() * out_radius_B;


                /* Conservative advancement, tighter bound using support functions. */
                double max_vel = clam::dot(shortest_dist_n, partB.vel);
                {
                    double abs_omega_A = partA.ang_vel.length();
                    if(abs_omega_A > 0.0){
                        auto inv_rot_A = partA.xform.rot_.inv();
                        clam::Vec3d c1  = inv_rot_A.rotate(clam::cross(partA.ang_vel, shortest_dist_n));
                        clam::Vec3d c1p = clam::cross(inv_rot_A.rotate(partA.ang_vel), c1) / abs_omega_A;
                        max_vel += clam::dot(c1, partA.xform.size_ * a.support(c1)) + clam::dot(c1p, partA.xform.size_ * a.support(c1p));
                    }
                    double abs_omega_B = partB.ang_vel.length();
                    if(abs_omega_B > 0.0){
                        auto inv_rot_B = partB.xform.rot_.inv();
                        clam::Vec3d c2  = inv_rot_B.rotate(clam::cross(partB.ang_vel, shortest_dist_n));
                        clam::Vec3d c2p = clam::cross(inv_rot_B.rotate(partB.ang_vel), c2) / abs_omega_B;
                        max_vel += clam::dot(c2, partB.xform.size_ * b.support(c2)) + clam::dot(c2p, partB.xform.size_ * b.support(c2p));
                    }
                }

                if(max_vel < 0.0) break;

                double max_advance = (distance - sim_.closest_distance_tol_) / max_vel;
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

                prev_distance = distance;
            }
        }

        return ParticleEvent::None();
    }

    ParticleEvent operator()(const shape::Sphere& a, const shape::Sphere& b)const{
        const Particle& partA = sim_.particles_[pa_idx_];
        const Particle& partB = sim_.particles_[pb_idx_];

        clam::Vec3d dist = sim_.pbc_.minImage(partB.xform.pos_ + partB.vel * (sim_.time_ - partB.time) - partA.xform.pos_);
        clam::Vec3d rel_vel = partA.vel - partB.vel;

        double time(0.0);
        clam::Vec3d normal;
        if(overlap::sphere_raycast(partA.xform.size_ * a.radius() + partB.xform.size_ * b.radius() + sim_.closest_distance_tol_, dist, rel_vel, time, &normal)){
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

        Transform bbox = sim_.boxes_[pid_];
        auto bbox_inv_rot = bbox.rot_.inv();

        Particle part = sim_.particles_[pid_];
        stream_rotation(part, sim_.time_);
        part.xform.pos_ = bbox_inv_rot.rotate(sim_.pbc_.minImage(part.xform.pos_ - bbox.pos_));
        part.xform.rot_ = bbox_inv_rot * part.xform.rot_;
        part.vel        = bbox_inv_rot.rotate(part.vel);
        part.ang_vel    = bbox_inv_rot.rotate(part.ang_vel);
        part.time       = 0.0;

        clam::Vec3d half_size = sim_.boxes_[pid_].size_ * boost::get<shape::Box>(*sim_.box_shapes_[part.shape_id]).extent() + sim_.obb_margin_;
        double out_radius = part.xform.size_ * shape_outradius(shape);

        while(true){
            double max_advance = std::numeric_limits<double>::max();
            double distance = 0.0;
            int max_id = 0;
            auto inv_rot = part.xform.rot_.inv();
            for(int i = 0; i < 3; ++i){
                double max_vel = part.vel[i] + part.ang_vel.length() * out_radius;
                if(max_vel > 0.0){
                    auto dir = clam::Vec3d(0.0);
                    dir[i] = 1.0;
                    auto dir_p = inv_rot.rotate(dir);
                    double dist = half_size[i] - part.xform.size_ * part.xform.rot_.rotate(shape.support(dir_p))[i] - part.xform.pos_[i];
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
                    double dist = half_size[i] + part.xform.size_ * part.xform.rot_.rotate(shape.support(dir_p))[i] + part.xform.pos_[i];
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
                    auto dir_p = part.xform.rot_.inv().rotate(dir);
                    distance = half_size[idx] - dir[idx] * (part.xform.size_ * part.xform.rot_.rotate(shape.support(dir_p))[idx] + part.xform.pos_[idx]);
                }

                return ParticleEvent::NeighborhoodCross(sim_.time_ + time, pid_, 0);
            }

            time += max_advance;

            stream_position(part, time);
            stream_rotation(part, time);
            part.time = time;
        }

        return ParticleEvent::None();
    }

    ParticleEvent operator()(const shape::Sphere& a)const{
        Particle part = sim_.particles_[pid_];

        clam::Vec3d rel_vel = -part.vel;
        clam::Vec3d dist    = sim_.pbc_.minImage(part.xform.pos_ - sim_.boxes_[pid_].pos_);

        double time(0.0);
        bool found = overlap::sphere_raycast_full(sim_.obb_margin_, dist, rel_vel, time);
        assert(found);
        _UNUSED(found);

        return ParticleEvent::NeighborhoodCross(sim_.time_ + time, pid_, 0);
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

ParticleEvent Simulation::get_neighborhood_cross_event(int pid)const{
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

    event_mgr_.push(pA, get_neighborhood_cross_event(pA));
    event_mgr_.push(pB, get_neighborhood_cross_event(pB));

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

void Simulation::run_neighborhood_cross_event(const ParticleEvent& event){
    int pid = event.pid_;

    for(int nb: nnl_[pid]){
        nnl_[nb].erase(std::find(nnl_[nb].begin(), nnl_[nb].end(), pid));
    }

    update_particle(particles_[pid], time_, pbc_);
    boxes_[pid] = particles_[pid].xform;

    cll_.update(pid, boxes_[pid].pos_);

    auto shape_id_a = particles_[pid].shape_id;

    auto bba = boxes_[pid];
    auto inv_rot = bba.rot_.inv();

    std::vector<size_t> nnl_new;

    for(int cid: cll_.cell_nbs(cll_.cell_index(pid))){
        for(int n: cll_.cell_content(cid)){
            if(n == pid) continue;
            auto shape_id_b = particles_[n].shape_id;
            auto bbb = boxes_[n];

            bbb.pos_ = inv_rot.rotate(pbc_.minImage(bbb.pos_ - bba.pos_));
            bbb.rot_ = inv_rot * bbb.rot_;

            if(overlap::bv_overlap(bba, *box_shapes_[shape_id_a], bbb, *box_shapes_[shape_id_b], obb_margin_)){
                nnl_new.push_back(n);
                nnl_[n].push_back(pid);
                if(std::find(nnl_[pid].begin(), nnl_[pid].end(), n) == nnl_[pid].end()){
                    auto event = get_collision_event(pid, n);
                    if(event.get_type() != PE_NONE) event_mgr_.push(pid, event);
                }
            }
        }
    }
    nnl_[pid] = nnl_new;

    event_mgr_.push(pid, get_neighborhood_cross_event(pid));
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

Simulation::~Simulation(void){
    for(size_t i = 0; i < shapes_.size(); ++i) delete box_shapes_[i];
    delete[] box_shapes_;
    delete[] boxes_;
    delete[] nnl_;
}

//TODO: Perhaps add exceptions to constructor
Simulation::Simulation(const Configuration& config):
    time_(0.0),
    statistics_start_time_(0.0),
    closest_distance_tol_(1.0e-6), //@note: increase tolerance to increase performance.
    obb_margin_(0.1),
    max_collision_time_(5.0),
    av_momentum_transfer_(0.0),
    av_kinetic_delta_(0.0), base_kinetic_energy_(0.0), kinetic_delta_(0.0),
    max_inflight_time_(0.0),
    n_collision_events_(0),
    config_(config),
    pbc_(config_.pbc_), particles_(config_.particles_), shapes_(config_.shapes_)
{
    mtGen_.seed(std::time(NULL));

    auto n_part = particles_.size();

    if(n_part) event_mgr_.resize(n_part);
    else return;

    //Initialize number of collisions to zero
    n_collisions_.resize(n_part, 0);

    //Initialize obb shapes
    box_shapes_ = new bounding_volume::Variant*[shapes_.size()];
    for(size_t i = 0; i < shapes_.size(); ++i) box_shapes_[i] = new bounding_volume::Variant(calc_bounding_volume(*shapes_[i]));

    boxes_ = new Transform[n_part];
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

        particle.xform.pos_ = pbc_.apply(particle.xform.pos_ + 0.1); //Make sure pbc are correct at start.
        boxes_[i] = particle.xform;

        //NOTE: Can cache radii for each shape
        double outradius = particle.xform.size_ * bv_outradius(*box_shapes_[particle.shape_id]) + obb_margin_;
        max_radius = std::max(max_radius, outradius);
    }
    sys_vel = sys_vel * (1.0 / n_part);

    double min_cell_size = 2.0 * max_radius + 0.01;

    //The box needs to be at least to cell sizes on each dimension
    //or else we might get overlaps, since we're not calculating
    //the collision time with all periodic images.
    assert(2.0 * min_cell_size < pbc_.getSize()[0]);
    assert(2.0 * min_cell_size < pbc_.getSize()[1]);
    assert(2.0 * min_cell_size < pbc_.getSize()[2]);

    //Initialize cell list
    cll_.init(n_part, pbc_.getSize(), min_cell_size);
    for(size_t i = 0; i < n_part; ++i){
        particles_[i].vel -= sys_vel;
        base_kinetic_energy_ += 0.5 * particles_[i].vel.length2();
        cll_.add(i, particles_[i].xform.pos_);
    }

    //Initialize nnl
    nnl_ = new std::vector<size_t>[n_part];

    foreach_pair(cll_, [this](int i, int j) -> bool {
        auto shape_id_i = particles_[i].shape_id;
        auto shape_id_j = particles_[j].shape_id;
        auto bbi = boxes_[i];
        auto bbj = boxes_[j];

        auto inv_rot = bbi.rot_.inv();
        bbj.pos_ = inv_rot.rotate(pbc_.minImage(boxes_[j].pos_ - boxes_[i].pos_));
        bbj.rot_ = inv_rot * bbj.rot_;

        if(overlap::bv_overlap(bbi, *box_shapes_[shape_id_i], bbj, *box_shapes_[shape_id_j], obb_margin_)){
            nnl_[i].push_back(j);
            nnl_[j].push_back(i);
            auto event = get_collision_event(i, j);
            if(event.get_type() != PE_NONE) event_mgr_.push(i, event);
        }
#ifndef NDEBUG
//Check for overlaps at start of simulation.
        Transform ta = particles_[i].xform;
        Transform tb = particles_[j].xform;
        tb.pos_ = pbc_.minImage(tb.pos_ - ta.pos_);
        ta.pos_ = 0.0;
        assert(overlap::shape_overlap(ta, *shapes_[particles_[i].shape_id], tb, *shapes_[particles_[j].shape_id]) == false);
#endif
        return false;
    });

    for(size_t i = 0; i < n_part; ++i){
        auto event = get_neighborhood_cross_event(i);
        assert(event.get_type() != PE_NONE);
        event_mgr_.push(i, event);
    }

    event_mgr_.init();
}

void Simulation::run(double end_time, PeriodicCallback& output_condition){

    is_running_ = true;

    while(is_running_){
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
        case PE_NEIGHBORHOOD_CROSS:
            run_neighborhood_cross_event(next_event);
            break;
        default:
            is_running_ = false;
            break;
        }

        output_condition(time_);

        if(time_ >= end_time) is_running_ = false;

        //Optimization: After we have collected enough events, we set the
        //maximum time for which to search for collisions to 1.5 times
        //the maximum in-flight time of each particle. Of course, it
        //might still be possible to miss a collision, but unlikely.
        if(n_collision_events_ == 40000) max_collision_time_ = 1.5 * max_inflight_time_;
    }
}

void Simulation::stop(void){
    is_running_ = false;
}

void Simulation::restart(void){
    reset_statistics();
    event_mgr_.clear();

    //The adress of nnl_ will be different if we deserialized earlier
    //TODO: We need a more reliable way to reseed, and if possible, it should
    //always be the same if rerun with the same parameters. Either that, or
    //it should be the user's responsiblity to reseed.
    mtGen_.seed(std::time(NULL) + reinterpret_cast<unsigned long>(nnl_));

    std::uniform_real_distribution<double> dist(-1.0, 1.0);
    clam::Vec3d sys_vel;
    for(size_t pid = 0; pid < particles_.size(); ++pid){
        nnl_[pid].clear();
        update_particle(particles_[pid], time_, pbc_);
        boxes_[pid] = particles_[pid].xform;
        cll_.update(pid, boxes_[pid].pos_);

        double x1, x2, r;
        do{
            x1 = dist(mtGen_);
            x2 = dist(mtGen_);
            r = x1 * x1 + x2 * x2;
        }while(r >= 1.0);
        double s = 2.0 * sqrt(1.0 - r);

        clam::Vec3d vec(x1 * s, x2 * s, 1.0 - 2.0 * r);
        sys_vel += vec;
        particles_[pid].vel = vec;
    }
    sys_vel = sys_vel * (1.0 / particles_.size());
    base_kinetic_energy_ = 0.0;
    for(size_t pid = 0; pid < particles_.size(); ++pid){
        particles_[pid].vel -= sys_vel;
        particles_[pid].ang_vel = 0.0;
        base_kinetic_energy_ += 0.5 * particles_[pid].vel.length2();
    }

    foreach_pair(cll_, [this](int i, int j) -> bool {
        auto shape_id_i = particles_[i].shape_id;
        auto shape_id_j = particles_[j].shape_id;
        auto bbi = boxes_[i];
        auto bbj = boxes_[j];

        auto inv_rot = bbi.rot_.inv();
        bbj.pos_ = inv_rot.rotate(pbc_.minImage(boxes_[j].pos_ - boxes_[i].pos_));
        bbj.rot_ = inv_rot * bbj.rot_;

        if(overlap::bv_overlap(bbi, *box_shapes_[shape_id_i], bbj, *box_shapes_[shape_id_j], obb_margin_)){
            nnl_[i].push_back(j);
            nnl_[j].push_back(i);
            auto event = get_collision_event(i, j);
            if(event.get_type() != PE_NONE) event_mgr_.push(i, event);
        }
#ifndef NDEBUG
        //There should be no overlaps!
        Transform ta = particles_[i].xform;
        Transform tb = particles_[j].xform;
        tb.pos_ = pbc_.minImage(tb.pos_ - ta.pos_);
        ta.pos_ = 0.0;
        assert(overlap::shape_overlap(ta, *shapes_[particles_[i].shape_id], tb, *shapes_[particles_[j].shape_id]) == false);
#endif
        return false;
    });

    for(size_t i = 0; i < particles_.size(); ++i){
        auto event = get_neighborhood_cross_event(i);
        assert(event.get_type() != PE_NONE);
        event_mgr_.push(i, event);
        event_mgr_.update(i);
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

bool Simulation::check_overlaps(void)const{
    for(size_t i = 0; i < particles_.size(); ++i){
        for(size_t j: nnl_[i]){
            if(j > i){
                Particle pa = particles_[i];
                Particle pb = particles_[j];
                update_particle(pa, time_, pbc_);
                update_particle(pb, time_, pbc_);
                pb.xform.pos_ = pbc_.minImage(pb.xform.pos_ - pa.xform.pos_);
                pa.xform.pos_ = 0.0;

                if(overlap::shape_overlap(pa.xform, *shapes_[particles_[i].shape_id], pb.xform, *shapes_[particles_[j].shape_id])){
                    return true;
                }
            }
        }
    }

    return false;
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
    serialize(ar, sim.obb_margin_);
    serialize(ar, sim.max_collision_time_);
    serialize(ar, sim.av_momentum_transfer_);
    serialize(ar, sim.av_kinetic_delta_);
    serialize(ar, sim.base_kinetic_energy_);
    serialize(ar, sim.kinetic_delta_);
    serialize(ar, sim.max_inflight_time_);
    serialize(ar, sim.n_collision_events_);
    serialize(ar, sim.n_collisions_);
    serialize(ar, sim.config_);

    for(size_t i = 0; i < sim.shapes_.size(); ++i) serialize(ar, *sim.box_shapes_[i]);
    serialize(ar, sim.boxes_, sim.particles_.size());
    for(size_t i = 0; i < sim.particles_.size(); ++i) serialize(ar, sim.nnl_[i]);

    serialize(ar, sim.event_mgr_);
    serialize(ar, sim.cll_);
}

void deserialize(Archive& ar, Simulation* sim){
    deserialize(ar, &sim->time_);
    deserialize(ar, &sim->prev_time_);
    deserialize(ar, &sim->statistics_start_time_);
    deserialize(ar, &sim->closest_distance_tol_);
    deserialize(ar, &sim->obb_margin_);
    deserialize(ar, &sim->max_collision_time_);
    deserialize(ar, &sim->av_momentum_transfer_);
    deserialize(ar, &sim->av_kinetic_delta_);
    deserialize(ar, &sim->base_kinetic_energy_);
    deserialize(ar, &sim->kinetic_delta_);
    deserialize(ar, &sim->max_inflight_time_);
    deserialize(ar, &sim->n_collision_events_);
    deserialize(ar, &sim->n_collisions_);
    deserialize(ar, &sim->config_);

    sim->box_shapes_ = new bounding_volume::Variant*[sim->shapes_.size()];
    for(size_t i = 0; i < sim->shapes_.size(); ++i) deserialize(ar, &sim->box_shapes_[i]);
    sim->boxes_ = new Transform[sim->particles_.size()];
    deserialize(ar, sim->boxes_, sim->particles_.size());

    sim->nnl_ = new std::vector<size_t>[sim->particles_.size()];
    for(size_t i = 0; i < sim->particles_.size(); ++i) deserialize(ar, &sim->nnl_[i]);

    deserialize(ar, &sim->event_mgr_);
    deserialize(ar, &sim->cll_);
}
