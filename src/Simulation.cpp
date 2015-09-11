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
        //double max_dist = max_overlap_distance(partA.size, *sim_.shapes_[partA.shape_id], partB.size, *sim_.shapes_[partB.shape_id]);
        //if(dist.length() <= max_dist){
        //    double time(0.0);
        //    clam::Vec3d relVel = partA.vel - partB.vel;
        //    if(overlap::sphere_raycast(max_dist, dist, relVel, time)){
        //        //PE_POSSIBLE_COLLISION
        //        return ParticleEvent(sim_.time_ + time, pa_idx_, sim_.n_part_ + 1 + pb_idx_, sim_.nCollisions_[pb_idx_]);
        //    }
        //}
        //else{
            partA.pos = 0.0;
            partB.pos = dist;
            clam::Vec3d relVel = partA.vel - partB.vel;
            double ispeed = 1.0 / relVel.length();
            double time = 10000.0;
            if(overlap::gjk_raycast(partA, a, partB, b, relVel * ispeed, time)){
                if(time > 0.0) return ParticleEvent(sim_.time_ + time * ispeed, pa_idx_, pb_idx_ + 1, sim_.nCollisions_[pb_idx_]);
            }
        //}

        return ParticleEvent();
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
    if(overlap::sphere_raycast(partA.size * a.radius() + partB.size * b.radius(), dist, relVel, time)){
        return ParticleEvent(sim_.time_ + time, pa_idx_, pb_idx_ + 1, sim_.nCollisions_[pb_idx_]);
    }
    else return ParticleEvent();
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
    return ParticleEvent(time + time_, pid, cellOffset + 2 * n_part_ + 1);
}

void Simulation::updateParticle(int pid){
    if(particles_[pid].time < time_){
        particles_[pid].pos += particles_[pid].vel * (time_ - particles_[pid].time);
        particles_[pid].pos  = pbc_.apply(particles_[pid].pos);
        particles_[pid].time = time_;
    }
}

//NOTE: For simplicity, for now we assume equal mass spheres
void Simulation::runCollisionEvent(const ParticleEvent& event){
    int pA = event.pid_;
    int pB = event.id_ - 1;

    if(nCollisions_[pB] != event.optional_){
        eventManager_.update(pA);
        return;
    }

    updateParticle(pA);
    updateParticle(pB);

    {
        clam::Vec3d relVel   = particles_[pA].vel - particles_[pB].vel;
        clam::Vec3d relPos   = pbc_.minImage(particles_[pA].pos - particles_[pB].pos);
        clam::Vec3d deltaVel = relPos * (clam::dot(relPos, relVel) / clam::dot(relPos, relPos));

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
                if(event.id_ != PE_NONE) eventManager_.push(pA, event);
            }
        }
    }
    for(int cid: cll_.getNeighbourIterator(pB)){
        for(int n: cll_.getCellIterator(cid)){
            if(n != pA && n != pB){
                //updateParticle(n);
                auto event = getCollisionEvent(pB, n);
                if(event.id_ != PE_NONE) eventManager_.push(pB, event);
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
    int pB = event.id_  - n_part_ - 1;

    if(nCollisions_[pB] != event.optional_){
        eventManager_.update(pA);
        return;
    }

    //should we?
    //updateParticle(pA);
    //updateParticle(pB);

    {
        Particle partA = particles_[pA];
        Particle partB = particles_[pB];

        partA.pos = 0.0;
        partB.pos = pbc_.minImage(partB.pos + partB.vel * (time_ - partB.time) - partA.pos);
        clam::Vec3d relVel = partA.vel - partB.vel;
        double ispeed = 1.0 / relVel.length();
        double time = 10000.0;
        if(overlap::shape_raycast(partA, *shapes_[partA.shape_id], partB, *shapes_[partA.shape_id], relVel * ispeed, time)){
            ParticleEvent event = ParticleEvent(time_ + time * ispeed, pA, pB + 1, nCollisions_[pB]);
            if(event.id_ != PE_NONE) eventManager_.push(pA, event);
        }
    }
}

void Simulation::runCellCrossEvent(const ParticleEvent& event){
    int pid     = event.pid_;
    int coffset = event.id_ - 2 * n_part_ - 1;
    cll_.move(pid, coffset);
    updateParticle(pid);
    for(int cid: cll_.getDirNeighbourIterator(pid, coffset)){
        for(int n: cll_.getCellIterator(cid)){
            //updateParticle(n);
            auto event = getCollisionEvent(pid, n);
            if(event.id_ != PE_NONE) eventManager_.push(pid, event);
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

    //Find initial collision events
    for(int i = 0; i < n_part_; ++i){
        for(int cid: cll_.getNeighbourIterator(i)){
            for(int j: cll_.getCellIterator(cid)){
                if(i != j){
                    auto event = getCollisionEvent(i, j);
                    if(event.id_ != PE_NONE) eventManager_.push(i, event);
                }
            }
        }
        auto event = getCellCrossEvent(i);
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
        outputCondition(nextEvent.time_);
        time_ = nextEvent.time_;
        //printf("%.16lf\n", time_);

        switch(nextEvent.getType(n_part_)){
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

