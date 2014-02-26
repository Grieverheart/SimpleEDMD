#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include "include/Simulation.h"
#include "include/EventManager.h"

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
            positions_.push_back(coords);
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
        fprintf(fp, "%f\t%f\t%f\t", positions_[i].x, positions_[i].y, positions_[i].z);
        fprintf(fp, "%f\n", radii_[i]);
    }
    fclose(fp);
}

void Simulation::addSphere(Vec3d pos, double radius){
    positions_.push_back(pos);
    radii_.push_back(radius);
    ++nSpheres_;
}

Vec3d Simulation::applyPeriodicBC(const Vec3d& vec)const{
    Vec3d retVec(0.0);
    for(size_t i = 0; i < 3; ++i) retVec[i] = remainder(vec[i], boxSize_);
    return retVec;
}

//WARNING: Needs re-evaluation
bool Simulation::raySphereIntersection(double radius, const Vec3d& pos, const Vec3d& dir, double& t)const{
    double dirInvLength = 1.0 / sqrt(dot(dir, dir));
    Vec3d dn  = dir * dirInvLength; //Normalize
    double s  = dot(pos, dn);
    double l2 = dot(pos, pos);
    double r2 = radius * radius;
    if(s < 0.0 && l2 > r2) return false;

    double m2 = l2 - s * s;
    if(m2 > r2) return false;

    double q = sqrt(r2 - m2);
    if(l2 > r2) t = s - q;
    else{
        t = s + q;
        printf("%f, %f\n", l2, r2);
    }
    t *= dirInvLength;

    return true;
}

CollisionEvent* Simulation::getCollisionEvent(size_t pA, size_t pB)const{
    Vec3d posA = positions_[pA],  posB = positions_[pB];
    Vec3d velA = velocities_[pA], velB = velocities_[pB];

    Vec3d dist   = applyPeriodicBC(posB - posA);
    Vec3d relVel = velA - velB;

    Time time(0.0);
    bool isCollision = raySphereIntersection(radii_[pA] + radii_[pB], dist, relVel, time);

    if(isCollision) return new CollisionEvent(time + time_, pA, pB, nCollisions_[pB]);
    else return nullptr;
}

void Simulation::updateParticle(size_t pID){
    if(times_[pID] != time_){
        positions_[pID] += velocities_[pID] * (time_ - times_[pID]);
        times_[pID] = time_;
    }
}

//NOTE: For simplicity, for now we assume equal mass spheres
void Simulation::runCollisionEvent(const CollisionEvent& event){
    size_t pA = event.pA;
    size_t pB = event.pB;

    if(nCollisions_[pB] != event.nBCollisions){
        eventManager_->update(pA);
        return;
    }

    updateParticle(pA);
    updateParticle(pB);

    {
        Vec3d relVel   = velocities_[pA] - velocities_[pB];
        Vec3d relPos   = applyPeriodicBC(positions_[pA] - positions_[pB]);
        Vec3d deltaVel = relPos * (dot(relPos, relVel) / dot(relPos, relPos));

        velocities_[pA] -= deltaVel;
        velocities_[pB] += deltaVel;
    }

    ++nCollisions_[pA];
    ++nCollisions_[pB];

    eventManager_->clear(pA);
    eventManager_->clear(pB);

    //Recalculate collision events
    for(size_t n = 0; n < nSpheres_; ++n){
        if(n != pA && n != pB){
            updateParticle(n);
            auto eventA = getCollisionEvent(pA, n);
            auto eventB = getCollisionEvent(pB, n);
            if(eventA) eventManager_->push(pA, eventA);
            if(eventB) eventManager_->push(pB, eventB);
        }
    }
    eventManager_->update(pA);
    eventManager_->update(pB);
}

bool Simulation::init(void){
    if(nSpheres_) eventManager_ = new EventManager(nSpheres_, 5000.0, 250000);
    else return false;

    //Initialize particle times to zero
    times_.resize(nSpheres_, 0.0);

    //Initialize number of collisions to zero
    nCollisions_.resize(nSpheres_, 0);

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
        velocities_.push_back(vec);
    }

    //Find initial collision events
    int nEvents = 0;
    for(size_t i = 0; i < nSpheres_; ++i){
        for(size_t j = i + 1; j < nSpheres_; ++j){
            auto event = getCollisionEvent(i, j);
            if(event) eventManager_->push(i, event);
        }
        if(!eventManager_->empty(i)) ++nEvents;
        eventManager_->update(i);
    }
    printf("%d\n", nEvents);

    return true;
}

void Simulation::run(void){

    Time endTime = 10.0;
    
    bool running = true;
    int nEvents = 0;
    Time snapTime = 0.1;
    while(running){
        Event* nextEvent = eventManager_->getNextEvent();
        time_ = nextEvent->time_;
        //std::cout << time_ << std::endl;

        switch(nextEvent->getType()){
        case EVT_COLLISION:{
                auto collisionEvent = static_cast<CollisionEvent*>(nextEvent);
                runCollisionEvent(*collisionEvent);
                delete collisionEvent;
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

