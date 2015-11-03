#include <cstdio>
#include <cmath>
#include <limits>
#include "Simulation.h"
#include "shape/variant.h"
#include "overlap/ray_casting.h"
#include "overlap/gjk.h"
#include "overlap/overlap.h"

//static void print_particle(const Particle& p){
//    printf("pos: glm::vec3(%e, %e, %e)\n", p.pos[0], p.pos[1], p.pos[2]);
//    clam::Vec3d axis;
//    double angle;
//    p.rot.toAxisAngle(angle, axis);
//    printf("rot: %e, glm::vec3(%e, %e, %e)\n", angle, axis[0], axis[1], axis[2]);
//    printf("vel: glm::vec3(%e, %e, %e)\n", p.vel[0], p.vel[1], p.vel[2]);
//    printf("ang_vel: glm::vec3(%e, %e, %e)\n", p.ang_vel[0], p.ang_vel[1], p.ang_vel[2]);
//}

namespace{
    template<typename T>
    T sqr(const T& x){
        return x * x;
    }

    class ShapeOutRadiusVisitor: public boost::static_visitor<double>{
    public:
        double operator()(const shape::Polyhedron& poly)const{
            return poly.out_radius();
        }

        double operator()(const shape::Sphere& sph)const{
            return sph.radius();
        }
    };

    inline double shape_outradius(const shape::Variant& shape){
        return boost::apply_visitor(ShapeOutRadiusVisitor(), shape);
    }

    inline double max_overlap_distance(double size_a, const shape::Variant& sa, double size_b, const shape::Variant& sb){
        return size_a * boost::apply_visitor(ShapeOutRadiusVisitor(), sa) +
               size_b * boost::apply_visitor(ShapeOutRadiusVisitor(), sb);
    }

    class ShapeCollisionResolutionVisitor: public boost::static_visitor<bool> {
    public:
        ShapeCollisionResolutionVisitor(Particle& pa, Particle& pb, const CubicPBC& pbc):
            pa_(pa), pb_(pb), pbc_(pbc)
        {}

        template<typename T, typename U>
        bool operator()(const T& a, const U& b)const{
            clam::Vec3d relPos = pbc_.minImage(pb_.pos - pa_.pos);
            Particle temp_pa = pa_;
            Particle temp_pb = pb_;
            temp_pb.pos = relPos;
            temp_pa.pos = 0.0;

            clam::Vec3d point_on_a, point_on_b;
            clam::Vec3d dist_vec = overlap::gjk_closest_points(temp_pa, a, temp_pb, b, point_on_a, point_on_b);
            if(dist_vec.length2() == 0.0) return false;

            //@note: We're using the distance vector itself for better accuracy.
            clam::Vec3d normal = -dist_vec / dist_vec.length();

            point_on_b -= relPos;
            clam::Vec3d relVel = pb_.vel - pa_.vel + clam::cross(pb_.ang_vel, point_on_b) - clam::cross(pa_.ang_vel, point_on_a);
            clam::Vec3d tangent_a = clam::cross(point_on_a, normal);
            clam::Vec3d tangent_b = clam::cross(point_on_b, normal);

            double momentum_delta = 
                2.0 * clam::dot(normal, relVel) /
                (2.0 + tangent_a.length2()  + tangent_b.length2());

            pa_.vel += momentum_delta * normal;
            pb_.vel -= momentum_delta * normal;

            pa_.ang_vel += momentum_delta * tangent_a;
            pb_.ang_vel -= momentum_delta * tangent_b;

            return true;
        }

        bool operator()(const shape::Sphere& a, const shape::Sphere& b)const{
            clam::Vec3d relPos = pbc_.minImage(pa_.pos - pb_.pos);
            double dist2 = relPos.length2();
            if(dist2 < sqr(a.radius() + b.radius())) return false;

            clam::Vec3d relVel = pa_.vel - pb_.vel;
            clam::Vec3d deltaVel = relPos * clam::dot(relPos, relVel) / dist2;

            pa_.vel -= deltaVel;
            pb_.vel += deltaVel;

            return true;
        }

    private:
        Particle& pa_;
        Particle& pb_;
        const CubicPBC& pbc_;
    };
}

void stream_position(Particle& particle, double time){
    particle.pos += particle.vel * (time - particle.time);
}

void stream_rotation(Particle& particle, double time){
    double ang_vel_abs = particle.ang_vel.length();
    if(ang_vel_abs != 0.0){
        particle.rot = clam::fromAxisAngle(
            ang_vel_abs * (time - particle.time),
            particle.ang_vel / ang_vel_abs
        ) * particle.rot;
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
        //TODO: Fix tolerance
        if(dist.length() > max_dist + 2.0 * sim_.closest_distance_tol_){
            double time(0.0);
            clam::Vec3d relVel = partA.vel - partB.vel;
            if(overlap::sphere_raycast(max_dist + 1.5 * sim_.closest_distance_tol_, dist, relVel, time)){
                return ParticleEvent::PossibleCollision(sim_.time_ + time, pa_idx_, pb_idx_, sim_.nCollisions_[pb_idx_]);
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

            double out_radius_A = partA.size * shape_outradius(a);
            double out_radius_B = partB.size * shape_outradius(b);

            //TODO: Look into high iteration number!!!
            while(true){
                clam::Vec3d shortest_dist = overlap::gjk_distance(partA, a, partB, b, sim_.closest_distance_tol_);

                clam::Vec3d shortest_dist_n = shortest_dist / shortest_dist.length();
                double max_vel = clam::dot(shortest_dist_n, partB.vel) +
                                 partA.ang_vel.length() * out_radius_A + partB.ang_vel.length() * out_radius_B;

                if(max_vel < 0.0) break;

                if(shortest_dist.length() < sim_.closest_distance_tol_){
                    int iter = 0;
                    while(shortest_dist.length() < sim_.closest_distance_tol_){
                        time *= 0.999;
                        stream_position(partB, time);
                        stream_rotation(partB, time);
                        stream_rotation(partA, time);
                        partB.time = time;
                        partA.time = time;
                        shortest_dist = overlap::gjk_distance(partA, a, partB, b, sim_.closest_distance_tol_);
                        //This should happen for grazing collisions
                        if(iter > 10) return ParticleEvent::None();
                        ++iter;
                    }

                    return ParticleEvent::Collision(sim_.time_ + time, pa_idx_, pb_idx_, sim_.nCollisions_[pb_idx_], shortest_dist);
                }

                double max_advance = shortest_dist.length() / max_vel;
                time += max_advance;

                //No collision.
                //NOTE: Why time > 1.0?
                if(time < 0.0 || time > 10.0) break;

                stream_position(partB, time);
                stream_rotation(partB, time);
                stream_rotation(partA, time);
                partB.time = time;
                partA.time = time;
            }
        }

        return ParticleEvent::None();
    }

private:
    const Simulation& sim_;
    int pa_idx_;
    int pb_idx_;
};

template<>
inline ParticleEvent Simulation::ShapeCollisionEventVisitor::operator()(const shape::Sphere& a, const shape::Sphere& b)const{
    const Particle& partA = sim_.particles_[pa_idx_];
    const Particle& partB = sim_.particles_[pb_idx_];

    clam::Vec3d dist = sim_.pbc_.minImage(partB.pos + partB.vel * (sim_.time_ - partB.time) - partA.pos);
    clam::Vec3d relVel = partA.vel - partB.vel;

    double time(0.0);
    clam::Vec3d normal;
    if(overlap::sphere_raycast(partA.size * a.radius() + partB.size * b.radius(), dist, relVel, time, &normal)){
        return ParticleEvent::Collision(sim_.time_ + time, pa_idx_, pb_idx_, sim_.nCollisions_[pb_idx_], normal);
    }
    else return ParticleEvent::None();
}

ParticleEvent Simulation::getCollisionEvent(int pA, int pB)const{
    //TODO: Pass correct shape
    return boost::apply_visitor(ShapeCollisionEventVisitor(*this, pA, pB), *shapes_[0], *shapes_[0]);
}

ParticleEvent Simulation::getCellCrossEvent(int pid)const{
    int cidx = cll_.getIndex(pid);
    double time(0.0);
    clam::Vec3d rpos = pbc_.minImage(particles_[pid].pos - cll_.getCellOrigin(cidx));
    int cellOffset = overlap::cell_raycast(cll_.getCellSize(), rpos, particles_[pid].vel, time);
    assert(time > 0.0);
    return ParticleEvent::CellCross(time_ + time, pid, cellOffset);
}

//TODO: We want to make this more general, and not a member function. e.g. it should take the
//particle, time, and the boundary condition as arguments.
void Simulation::updateParticle(int pid){
    if(particles_[pid].time < time_){
        stream_position(particles_[pid], time_);
        stream_rotation(particles_[pid], time_);
        particles_[pid].pos = pbc_.apply(particles_[pid].pos);
        particles_[pid].time = time_;
    }
}

//NOTE: For simplicity, for now we assume equal mass spheres
void Simulation::runCollisionEvent(const ParticleEvent& event){
    int pA = event.pid_;
    int pB = event.get_id();

    if(nCollisions_[pB] != event.optional_){
        eventManager_.update(pA);
        return;
    }

    updateParticle(pA);
    updateParticle(pB);

    //Resolve collision.
    bool resolved = boost::apply_visitor(
        ShapeCollisionResolutionVisitor(particles_[pA], particles_[pB], pbc_),
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
        auto event = getCollisionEvent(pA, pB);
        if(event.get_type() != PE_NONE){
            eventManager_.push(pA, event);
            eventManager_.update(pA);
        }
        return;
    }

    ++nCollisions_[pA];
    ++nCollisions_[pB];

    eventManager_.clear(pA);
    eventManager_.clear(pB);

    {
        auto event = getCollisionEvent(pA, pB);
        if(event.get_type() != PE_NONE){
            if(event.time_ > time_) eventManager_.push(pA, event);
            else assert(false);
        }
    }

    //Recalculate collision events
    for(int cid: cll_.getNeighbourIterator(pA)){
        for(int n: cll_.getCellIterator(cid)){
            if(n != pA && n != pB){
                //updateParticle(n);
                auto event = getCollisionEvent(pA, n);
                if(event.get_type() != PE_NONE) eventManager_.push(pA, event);
            }
        }
    }

    for(int cid: cll_.getNeighbourIterator(pB)){
        for(int n: cll_.getCellIterator(cid)){
            if(n != pA && n != pB){
                //updateParticle(n);
                auto event = getCollisionEvent(pB, n);
                if(event.get_type() != PE_NONE) eventManager_.push(pB, event);
            }
        }
    }

    {
        auto event = getCellCrossEvent(pA);
        assert(event.get_type() != PE_NONE);
        assert(event.time_ > 0.0);
        eventManager_.push(pA, event);
    }
    {
        auto event = getCellCrossEvent(pB);
        assert(event.get_type() != PE_NONE);
        assert(event.time_ > 0.0);
        eventManager_.push(pB, event);
    }

    eventManager_.update(pA);
    eventManager_.update(pB);
}

void Simulation::runPossibleCollisionEvent(const ParticleEvent& event){
    int pA = event.pid_;
    int pB = event.get_id();

    updateParticle(pA);

    if(nCollisions_[pB] != event.optional_){
        eventManager_.update(pA);
        return;
    }

    auto new_event = getCollisionEvent(pA, pB);
    if(new_event.get_type() != PE_NONE){
        assert(new_event.time_ > time_);
        assert(new_event.get_type() != PE_POSSIBLE_COLLISION);
        eventManager_.push(pA, new_event);
        eventManager_.update(pA);
    }
}

void Simulation::runCellCrossEvent(const ParticleEvent& event){
    int pid     = event.pid_;
    int coffset = event.get_id();
    cll_.move(pid, coffset);
    updateParticle(pid);
    for(int cid: cll_.getDirNeighbourIterator(pid, coffset)){
        for(int n: cll_.getCellIterator(cid)){
            //updateParticle(n);
            auto event = getCollisionEvent(pid, n);
            if(event.get_type() != PE_NONE) eventManager_.push(pid, event);
        }
    }
    eventManager_.push(pid, getCellCrossEvent(pid));
    eventManager_.update(pid);
}

bool Simulation::init(void){
    if(n_part_) eventManager_.resize(n_part_);
    else return false;

    //Initialize number of collisions to zero
    nCollisions_.resize(n_part_, 0);

    //Initialize cell list
    double max_radius = 0.0;
    for(auto particle: particles_){
        double outradius = particle.size * boost::apply_visitor(ShapeOutRadiusVisitor(), *shapes_[particle.shape_id]);
        max_radius = std::max(max_radius, outradius);
    }
    cll_.init(n_part_, pbc_.getSize(), 2.0 * max_radius + 0.01);

    //Initialize paricle velocities
    std::uniform_real_distribution<double> dist(-1.0, 1.0);
    for(int i = 0; i < n_part_; ++i){
        double x1, x2, r;
        do{
            x1 = dist(mtGen_);
            x2 = dist(mtGen_);
            r = x1 * x1 + x2 * x2;
        }while(r >= 1.0);
        double s = 2.0 * sqrt(1.0 - r);
        clam::Vec3d vec(x1 * s, x2 * s, 1.0 - 2.0 * r);

        systemVelocity_  += vec;
        particles_[i].vel = vec;
        cll_.add(i, particles_[i].pos);
    }
    systemVelocity_ = systemVelocity_ * (1.0 / n_part_);

    for(int i = 0; i < n_part_; ++i){
        for(int j = i + 1; j < n_part_; ++j){
            auto partA = particles_[i];
            auto partB = particles_[j];
            partB.pos = pbc_.minImage(partB.pos + partB.vel * (time_ - partB.time) - partA.pos);
            partA.pos = 0.0;
        }
    }

    //Find initial collision events
    //NOTE: This is wrong, we need through each pair only once!!!
    //for(int i = 0; i < n_part_; ++i){
    //    for(int cid: cll_.getNeighbourIterator(i)){
    //        for(int j: cll_.getCellIterator(cid)){
    //            if(i != j){
    //                auto event = getCollisionEvent(i, j);
    //                if(event.get_id() != PE_NONE) eventManager_.push(i, event);
    //            }
    //        }
    //    }
    //    auto event = getCellCrossEvent(i);
    //    eventManager_.push(i, event);
    //}
    //TODO: Use the cell list properly.
    for(int i = 0; i < n_part_; ++i){
        for(int j = i + 1; j < n_part_; ++j){
            auto event = getCollisionEvent(i, j);
            if(event.get_id() != PE_NONE){
                assert(event.time_ > 0.0);
                eventManager_.push(i, event);
            }
        }
        auto event = getCellCrossEvent(i);
        assert(event.time_ > 0.0);
        eventManager_.push(i, event);
    }
    eventManager_.init();

    return true;
}

void Simulation::run(double endTime, PeriodicCallback& outputCondition){

    bool running = true;
    unsigned int nEvents = 0;
    while(running){
        ParticleEvent nextEvent = eventManager_.getNextEvent();
        prev_time_ = time_;
        time_ = nextEvent.time_;
        printf("%.16lf, %d, %d\n", time_, nEvents, nextEvent.get_type());

        switch(nextEvent.get_type()){
        case PE_COLLISION:
            ++nEvents;
            runCollisionEvent(nextEvent);
            break;
        case PE_POSSIBLE_COLLISION:
            runPossibleCollisionEvent(nextEvent);
            break;
        case PE_CELLCROSS:
            runCellCrossEvent(nextEvent);
            break;
        default:
            running = false;
            break;
        }

        outputCondition(nextEvent.time_);

        if(time_ >= endTime) running = false;
    }
}

