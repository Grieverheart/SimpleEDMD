#include <cstdio>
#include <cmath>
#include "include/Simulation.h"
#include "include/EventManager.h"
#include "include/ray_intersections.h"

void Simulation::saveConfig(const char* filename){
    FILE *fp = fopen(filename, "w");
    fprintf(fp, "%d\n%f\n", nSpheres_, pbc_.getSize());

    for(int i = 0; i < nSpheres_; ++i){
        updateParticle(i);
        fprintf(fp, "%f\t%f\t%f\t", particles_[i].pos.x, particles_[i].pos.y, particles_[i].pos.z);
        fprintf(fp, "%f\n", particles_[i].radius);
    }
    fclose(fp);
}

ParticleEvent Simulation::getCollisionEvent(int pA, int pB)const{
    const Particle& partA = particles_[pA];
    const Particle& partB = particles_[pB];

    Vec3d dist   = pbc_.minImage(partB.pos + partB.vel * (time_ - partB.time) - partA.pos);
    Vec3d relVel = partA.vel - partB.vel;

    double time(0.0);
    bool isCollision = raySphereIntersection(partA.radius + partB.radius, dist, relVel, time);

    if(isCollision) return ParticleEvent(time + time_, pA, pB + 1, nCollisions_[pB]);
    else return ParticleEvent();
}

ParticleEvent Simulation::getCellCrossEvent(int pid)const{
    int cidx = cll_.getIndex(pid);
    double time(0.0);
    Vec3d rpos = pbc_.minImage(particles_[pid].pos - cll_.getCellOrigin(cidx));
    int cellOffset = rayCellIntersection(cll_.getCellSize(), rpos, particles_[pid].vel, time);
    return ParticleEvent(time + time_, pid, cellOffset + nSpheres_ + 1);
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
        Vec3d relVel   = particles_[pA].vel - particles_[pB].vel;
        Vec3d relPos   = pbc_.minImage(particles_[pA].pos - particles_[pB].pos);
        Vec3d deltaVel = relPos * (dot(relPos, relVel) / dot(relPos, relPos));

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
    int coffset = event.id_ - nSpheres_ - 1;
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
    for(auto particle: particles_) max_radius = std::max(max_radius, particle.radius);
    cll_.init(nSpheres_, pbc_.getSize(), 2.0 * max_radius);

    //Initialize paricle velocities
    std::uniform_real_distribution<double> dist(-1.0, 1.0);
    for(int i = 0; i < nSpheres_; ++i){
        Vec3d vec;
        double sum = 0.0;
        for(int j = 0; j < 3; ++j){
            vec[j] = dist(mtGen_);
            sum += vec[j] * vec[j];
        }
        vec = vec * (1.0 / sqrt(sum));
        particles_[i].vel = vec;
        cll_.add(i, particles_[i].pos);
    }

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

void Simulation::run(void){

    double endTime = 100.0;
    
    bool running = true;
    int nEvents = 0;
    double snapTime = 100.0;
    while(running){
        ParticleEvent nextEvent = eventManager_.getNextEvent();
        time_ = nextEvent.time_;
        //std::cout << time_ << std::endl;

        switch(nextEvent.getType(nSpheres_)){
        case PE_COLLISION:
            runCollisionEvent(nextEvent);
            break;
        case PE_CELLCROSS:
            runCellCrossEvent(nextEvent);
            break;
        default:
            running = false;
            break;
        }
        if(time_ >= snapTime){
            char buff[64];
            sprintf(buff, "Data/file%06d.dat", ++nEvents);
            saveConfig(buff);
            snapTime += 1.0;
        }
        if(time_ > endTime) break;
    }
}

