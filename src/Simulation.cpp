#include <cstdio>
#include <cmath>
#include <limits>
#include "Simulation.h"
#include "EventManager.h"
#include "shape/variant.h"
#include "overlap/ray_casting.h"
#include "overlap/gjk.h"
#include "overlap/overlap.h"

namespace{
    class ShapeOutRadiusVisitor: public boost::static_visitor<double>{
    public:
        double operator()(const shape::Polyhedron& poly)const{
            return poly.out_radius();
        }

        double operator()(const shape::Sphere& sph)const{
            return sph.radius();
        }
    };

    inline double max_overlap_distance(double size_a, const shape::Variant& sa, double size_b, const shape::Variant& sb){
        return size_a * boost::apply_visitor(ShapeOutRadiusVisitor(), sa) +
               size_b * boost::apply_visitor(ShapeOutRadiusVisitor(), sb);
    }
}

static inline double shape_outradius(const shape::Variant& shape){
    return boost::apply_visitor(ShapeOutRadiusVisitor(), shape);
}

void stream_position(Particle& particle, double time){
    particle.pos += particle.vel * (time - particle.time);
}

void stream_rotation(Particle& particle, double time){
    double ang_vel_abs = particle.ang_vel.length();
    if(ang_vel_abs > 0.0){
        particle.rot = clam::fromAxisAngle(
            ang_vel_abs * (time - particle.time),
            particle.ang_vel / ang_vel_abs
        ) * particle.rot;
    }
}

//static void print_particle(const Particle& p){
//    printf("pos: glm::vec3(%f, %f, %f)\n", p.pos[0], p.pos[1], p.pos[2]);
//    clam::Vec3d axis;
//    double angle;
//    p.rot.toAxisAngle(angle, axis);
//    printf("rot: %f, glm::vec3(%f, %f, %f)\n", angle, axis[0], axis[1], axis[2]);
//    printf("vel: glm::vec3(%f, %f, %f)\n", p.vel[0], p.vel[1], p.vel[2]);
//    printf("ang_vel: glm::vec3(%f, %f, %f)\n", p.ang_vel[0], p.ang_vel[1], p.ang_vel[2]);
//}

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
        if(dist.length() > max_dist){
            double time(0.0);
            clam::Vec3d relVel = partA.vel - partB.vel;
            if(overlap::sphere_raycast(max_dist, dist, relVel, time)){
                //PE_POSSIBLE_COLLISION
                return ParticleEvent::PossibleCollision(sim_.time_ + time, pa_idx_, pb_idx_, sim_.nCollisions_[pb_idx_]);
            }
        }
        else{
            double time(0.0);

            partB.pos  = dist;
            partB.vel  = partB.vel - partA.vel;
            partB.time = 0.0;
            partA.pos  = 0.0;
            partA.vel  = 0.0;
            partA.time = 0.0;

            double out_radius_A = partA.size * shape_outradius(a);
            double out_radius_B = partB.size * shape_outradius(b);

            while(true){
                clam::Vec3d shortest_dist = overlap::gjk_distance(partA, a, partB, b);
                clam::Vec3d shortest_dist_n = shortest_dist / shortest_dist.length();
                double max_vel = clam::dot(shortest_dist_n, partB.vel) +
                                 partA.ang_vel.length() * out_radius_A + partB.ang_vel.length() * out_radius_B;
                double max_advance = shortest_dist.length() / max_vel;
                max_advance -= 1.0e-8; //@note: We need to subtract a small amount to avoid overlapping.

                //No collision.
                if(max_advance < 0.0) break;

                if(shortest_dist.length2() < sim_.closest_distance_tol2_){
                    //@note: We could avoid calculating closest points and call closest_points
                    //once when we process the event. Do that first, and then we can check how
                    //to modify the event structure to accomondate the contact information.
                    //
                    //clam::Vec3d pa, pb;
                    //overlap::gjk_closest_points(partA, a, partB, b, pa, pb);
                    return ParticleEvent::Collision(sim_.time_ + time, pa_idx_, pb_idx_, sim_.nCollisions_[pb_idx_], shortest_dist);
                }
                time += max_advance;
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
    return ParticleEvent::CellCross(time + time_, pid, cellOffset);
}

//TODO: We want to make this more general, and not a member function. e.g. it should take the
//particle, time, and the boundary condition as arguments.
void Simulation::updateParticle(int pid){
    if(particles_[pid].time < time_){
        stream_position(particles_[pid], time_);
        stream_rotation(particles_[pid], time_);
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
    {
        //TODO: Here we need a function that gets us contact information.
        //It's nice to do it here because we save space on events, and
        //we can also immediately check if the event is valid due to 
        //precision problems.
        clam::Vec3d relVel   = particles_[pA].vel - particles_[pB].vel;
        clam::Vec3d deltaVel = event.normal_ * clam::dot(event.normal_, relVel);

        particles_[pA].vel -= deltaVel;
        particles_[pB].vel += deltaVel;
    }

    ++nCollisions_[pA];
    ++nCollisions_[pB];

    eventManager_.clear(pA);
    eventManager_.clear(pB);

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
    eventManager_.push(pA, getCellCrossEvent(pA));
    eventManager_.push(pB, getCellCrossEvent(pB));

    eventManager_.update(pA);
    eventManager_.update(pB);
}

void Simulation::runPossibleCollisionEvent(const ParticleEvent& event){
    int pA = event.pid_;
    int pB = event.get_id();


    if(nCollisions_[pB] != event.optional_){
        eventManager_.update(pA);
        return;
    }

    {
        Particle partA = particles_[pA];
        Particle partB = particles_[pB];

        partB.pos = pbc_.minImage(partB.pos + partB.vel * (time_ - partB.time) - partA.pos - partA.vel * (time_ - partA.time));
        partA.pos = 0.0;
        clam::Vec3d relVel = partA.vel - partB.vel;
        double ispeed = 1.0 / relVel.length();
        double time = 10000.0;
        clam::Vec3d normal;
        if(overlap::shape_raycast(partA, *shapes_[partA.shape_id], partB, *shapes_[partB.shape_id], relVel * ispeed, time, normal)){
            ParticleEvent event = ParticleEvent::Collision(time_ + time * ispeed, pA, pB, nCollisions_[pB], normal);
            eventManager_.push(pA, event);
        }
    }

    eventManager_.update(pA);
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
            if(event.get_id() != PE_NONE) eventManager_.push(pid, event);
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
            auto distance = overlap::shape_distance(partA, *shapes_[particles_[i].shape_id], partB, *shapes_[particles_[j].shape_id]);
            if(distance.length() <= 0.0) printf("%d, %d\n", i, j);
        }
    }

    //Find initial collision events
    for(int i = 0; i < n_part_; ++i){
        for(int cid: cll_.getNeighbourIterator(i)){
            for(int j: cll_.getCellIterator(cid)){
                if(i != j){
                    auto event = getCollisionEvent(i, j);
                    if(event.get_id() != PE_NONE) eventManager_.push(i, event);
                }
            }
        }
        auto event = getCellCrossEvent(i);
        eventManager_.push(i, event);
    }
    eventManager_.init();


    //auto partA = particles_[1279];
    //auto partB = particles_[855];
    //auto a = *shapes_[partA.shape_id];
    //auto b = *shapes_[partB.shape_id];

    //printf("partA: %d\n", 1279);
    //print_particle(partA);
    //printf("partB: %d\n", 855);
    //print_particle(partB);

    //clam::Vec3d shortest_dist = overlap::shape_distance(partA, a, partB, b);
    //printf("== %e ==\n", shortest_dist.length());

    return true;
}

void Simulation::run(double endTime, PeriodicCallback& outputCondition){

    bool running = true;
    unsigned int nEvents = 0;
    while(running){
        ParticleEvent nextEvent = eventManager_.getNextEvent();
        outputCondition(nextEvent.time_);
        time_ = nextEvent.time_;
        //printf("%.16lf\n", time_);

        switch(nextEvent.get_type()){
        case PE_COLLISION:
            ++nEvents;
            runCollisionEvent(nextEvent);
            break;
        case PE_POSSIBLE_COLLISION:
            ++nEvents;
            runPossibleCollisionEvent(nextEvent);
            break;
        case PE_CELLCROSS:
            runCellCrossEvent(nextEvent);
            break;
        default:
            running = false;
            break;
        }
        if(time_ >= endTime) running = false;
    }
}

