#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include "include/Simulation.h"
#include "include/EventManager.h"
#include "include/ray_intersections.h"

void Simulation::readConfig(const char* filename){
	char line[128];
    const char* delim = "\t";
	char *t1 = NULL;
	int i = 1,u = 0;
	
	FILE *fp;
	fp = fopen(filename,"r");
	if(!fp){
		printf("Couldn't open file %s\n",filename);
		return;
	}
	
	while(fgets(line,sizeof(line),fp) != NULL){
		u=0;
		if(line[strlen(line)-1] == '\n') line[strlen(line)-1] = '\0'; // Remove the niewline from the end of the string
        if(i == 1) nSpheres_ = atoi(line);
        else if(i == 2) boxSize_ = atof(line);
        else if(i > 2){
            Vec3d coords(0.0);
            for(t1 = strtok(line,delim); t1 != NULL; t1 = strtok(NULL, delim)){
                if(u<3)	coords[u] = atof(t1);
                else radii_.push_back(atof(t1));
                u++;
            }
            Particle part;
            part.pos = coords;
            particles_.emplace_back(part);
		}
		i++;
	}
	free(t1);
	fclose(fp);
}

void Simulation::saveConfig(const char* filename){
    FILE *fp = fopen(filename, "w");
    fprintf(fp, "%lu\n%f\n", nSpheres_, boxSize_);

    for(size_t i = 0; i < nSpheres_; ++i){
        updateParticle(i);
        fprintf(fp, "%f\t%f\t%f\t", particles_[i].pos.x, particles_[i].pos.y, particles_[i].pos.z);
        fprintf(fp, "%f\n", radii_[i]);
    }
    fclose(fp);
}

//Vec3d Simulation::applyPeriodicBC(const Vec3d& vec)const{
//    Vec3d retVec(0.0);
//    for(size_t i = 0; i < 3; ++i) retVec[i] = remainder(vec[i], boxSize_);
//    return retVec;
//}

//// Most efficient but vec is required to be -B < vec < B
//Vec3d Simulation::applyPeriodicBC(const Vec3d& vec)const{
//    Vec3d retVec(0.0);
//    for(size_t i = 0; i < 3; ++i){
//        retVec[i] = vec[i] - int(vec[i] * (2.0 / boxSize_)) * boxSize_;
//    }
//    return retVec;
//}

// More efficient but vec is required to be -1.5B < vec < 1.5B
Vec3d Simulation::applyPeriodicBC(const Vec3d& vec)const{
    Vec3d retVec(0.0);
    static double a = 2.0 / boxSize_;
    for(size_t i = 0; i < 3; ++i){
        int k = vec[i] * a;
        retVec[i] = vec[i] - k * boxSize_;
        k = retVec[i] * a;
        retVec[i] = retVec[i] - k * boxSize_;
    }
    return retVec;
}

CollisionEvent* Simulation::getCollisionEvent(size_t pA, size_t pB)const{
    const Particle& partA = particles_[pA];
    const Particle& partB = particles_[pB];

    Vec3d dist   = applyPeriodicBC(partB.pos - partA.pos);
    Vec3d relVel = partA.vel - partB.vel;

    double time(0.0);
    bool isCollision = raySphereIntersection(radii_[pA] + radii_[pB], dist, relVel, time);

    if(isCollision) return new CollisionEvent(time + time_, pA, pB, nCollisions_[pB]);
    else return nullptr;
}

CellCrossEvent* Simulation::getCellCrossEvent(size_t pid)const{
    int cidx = cll_.getIndex(pid);
    double time(0.0);
    Vec3d rpos = applyPeriodicBC(particles_[pid].pos - cll_.getCellOrigin(cidx));
    int cellOffset = rayCellIntersection(cll_.getCellSize(), rpos, particles_[pid].vel, time);
    return new CellCrossEvent(time + time_, pid, cellOffset);
}

void Simulation::updateParticle(size_t pid){
    if(particles_[pid].time != time_){
        particles_[pid].pos += particles_[pid].vel * (time_ - particles_[pid].time);
        for(int i = 0; i < 3; ++i) particles_[pid].pos[i] -= int(particles_[pid].pos[i] * (2.0 / boxSize_) - 1.0) * boxSize_;
        particles_[pid].time = time_;
    }
}

//NOTE: For simplicity, for now we assume equal mass spheres
void Simulation::runCollisionEvent(const CollisionEvent& event){
    size_t pA = event.pA;
    size_t pB = event.pB;

    if(nCollisions_[pB] != event.nBCollisions){
        eventManager_.update(pA);
        return;
    }

    updateParticle(pA);
    updateParticle(pB);

    {
        Vec3d relVel   = particles_[pA].vel - particles_[pB].vel;
        Vec3d relPos   = applyPeriodicBC(particles_[pA].pos - particles_[pB].pos);
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
        for(size_t n: cll_.getCellIterator(cid)){
            if(n != pA && n != pB){
                updateParticle(n);
                auto event = getCollisionEvent(pA, n);
                if(event) eventManager_.push(pA, event);
            }
        }
    }
    for(int cid: cll_.getNeighbourIterator(pB)){
        for(size_t n: cll_.getCellIterator(cid)){
            if(n != pA && n != pB){
                updateParticle(n);
                auto event = getCollisionEvent(pB, n);
                if(event) eventManager_.push(pB, event);
            }
        }
    }
    eventManager_.push(pA, getCellCrossEvent(pA));
    eventManager_.push(pB, getCellCrossEvent(pB));

    eventManager_.update(pA);
    eventManager_.update(pB);
}

void Simulation::runCellCrossEvent(const CellCrossEvent& event){
    size_t pid  = event.pid;
    int coffset = event.coffset;
    cll_.move(pid, coffset);
    updateParticle(pid);
    for(int cid: cll_.getDirNeighbourIterator(pid, coffset)){
        for(size_t n: cll_.getCellIterator(cid)){
            updateParticle(n);
            auto event = getCollisionEvent(pid, n);
            if(event) eventManager_.push(pid, event);
        }
    }
    eventManager_.push(pid, getCellCrossEvent(pid));
    eventManager_.update(pid);
}

bool Simulation::init(void){
    if(nSpheres_) eventManager_.resize(nSpheres_);
    else return false;

    particles_.resize(nSpheres_);

    //Initialize number of collisions to zero
    nCollisions_.resize(nSpheres_, 0);

    //Initialize cell list
    double max_radius = 0.0;
    for(auto radius: radii_) max_radius = std::max(max_radius, radius);
    cll_.init(nSpheres_, boxSize_, 2.0 * max_radius);

    //Initialize paricle velocities
    std::uniform_real_distribution<double> dist(-1.0, 1.0);
    for(size_t i = 0; i < nSpheres_; ++i){
        Vec3d vec(0.0);
        double sum = 0.0;
        for(size_t j = 0; j < 3; ++j){
            vec[j] = dist(mtGen_);
            sum += vec[j] * vec[j];
        }
        vec = vec * (1.0 / sqrt(sum));
        particles_[i].vel = vec;
        cll_.add(i, particles_[i].pos);
    }

    //Find initial collision events
    for(size_t i = 0; i < nSpheres_; ++i){
        for(int cid: cll_.getNeighbourIterator(i)){
            for(size_t j: cll_.getCellIterator(cid)){
                if(i != j){
                    auto event = getCollisionEvent(i, j);
                    if(event) eventManager_.push(i, event);
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

    double endTime = 10.0;
    
    bool running = true;
    int nEvents = 0;
    double snapTime = 0.1;
    while(running){
        Event* nextEvent = eventManager_.getNextEvent();
        time_ = nextEvent->time_;
        //std::cout << time_ << std::endl;

        switch(nextEvent->getType()){
        case EVT_COLLISION:{
                auto collisionEvent = static_cast<CollisionEvent*>(nextEvent);
                runCollisionEvent(*collisionEvent);
                delete collisionEvent;
            }
            break;
        case EVT_CELLCROSS:{
                auto cellEvent = static_cast<CellCrossEvent*>(nextEvent);
                runCellCrossEvent(*cellEvent);
                delete cellEvent;
            }
            break;
        default:
            running = false;
            break;
        }
        if(time_ > endTime) break;
        if(time_ >= snapTime){
            char buff[64];
            sprintf(buff, "Data/file%06d.dat", ++nEvents);
            saveConfig(buff);
            snapTime += 0.1;
        }
    }
}

