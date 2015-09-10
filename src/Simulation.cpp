#include <cstdio>
#include <cmath>
#include "Simulation.h"
#include "EventManager.h"
#include "ray_intersections.h"

//TODO: This will check for and return a PE_POSSIBLE_COLLISION in case the
//polyhedra are farther apart than the sum of their circumscribed radii.
ParticleEvent Simulation::getCollisionEvent(int pA, int pB)const{
    const Particle& partA = particles_[pA];
    const Particle& partB = particles_[pB];

    clam::Vec3d dist   = pbc_.minImage(partB.pos + partB.vel * (time_ - partB.time) - partA.pos);
    clam::Vec3d relVel = partA.vel - partB.vel;

    double time(0.0);
    bool isCollision = raySphereIntersection(partA.size + partB.size, dist, relVel, time);

    if(isCollision) return ParticleEvent(time + time_, pA, pB + 1, nCollisions_[pB]);
    else return ParticleEvent();
}

ParticleEvent Simulation::getCellCrossEvent(int pid)const{
    int cidx = cll_.getIndex(pid);
    double time(0.0);
    clam::Vec3d rpos = pbc_.minImage(particles_[pid].pos - cll_.getCellOrigin(cidx));
    int cellOffset = rayCellIntersection(cll_.getCellSize(), rpos, particles_[pid].vel, time);
    return ParticleEvent(time + time_, pid, cellOffset + 2 * nSpheres_ + 1);
}

void Simulation::updateParticle(int pid){
    if(particles_[pid].time != time_){
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
                if(event.id_ != 0) eventManager_.push(pA, event);
            }
        }
    }
    for(int cid: cll_.getNeighbourIterator(pB)){
        for(int n: cll_.getCellIterator(cid)){
            if(n != pA && n != pB){
                //updateParticle(n);
                auto event = getCollisionEvent(pB, n);
                if(event.id_ != 0) eventManager_.push(pB, event);
            }
        }
    }
    eventManager_.push(pA, getCellCrossEvent(pA));
    eventManager_.push(pB, getCellCrossEvent(pB));

    eventManager_.update(pA);
    eventManager_.update(pB);
}

void Simulation::runCellCrossEvent(const ParticleEvent& event){
    int pid     = event.pid_;
    int coffset = event.id_ - 2 * nSpheres_ - 1;
    cll_.move(pid, coffset);
    updateParticle(pid);
    for(int cid: cll_.getDirNeighbourIterator(pid, coffset)){
        for(int n: cll_.getCellIterator(cid)){
            //updateParticle(n);
            auto event = getCollisionEvent(pid, n);
            if(event.id_ != 0) eventManager_.push(pid, event);
        }
    }
    eventManager_.push(pid, getCellCrossEvent(pid));
    eventManager_.update(pid);
}

bool Simulation::init(void){
    if(nSpheres_) eventManager_.resize(nSpheres_);
    else return false;

    //Initialize number of collisions to zero
    nCollisions_.resize(nSpheres_, 0);

    //Initialize cell list
    double max_radius = 0.0;
    for(auto particle: particles_) max_radius = std::max(max_radius, particle.size);
    cll_.init(nSpheres_, pbc_.getSize(), 2.0 * max_radius);

    //Initialize paricle velocities
    std::uniform_real_distribution<double> dist(-1.0, 1.0);
    for(int i = 0; i < nSpheres_; ++i){
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
    systemVelocity_ = systemVelocity_ * (1.0 / nSpheres_);

    //Find initial collision events
    for(int i = 0; i < nSpheres_; ++i){
        for(int cid: cll_.getNeighbourIterator(i)){
            for(int j: cll_.getCellIterator(cid)){
                if(i != j){
                    auto event = getCollisionEvent(i, j);
                    if(event.id_ != 0) eventManager_.push(i, event);
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

        switch(nextEvent.getType(nSpheres_)){
        case PE_COLLISION:
            ++nEvents;
            runCollisionEvent(nextEvent);
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

