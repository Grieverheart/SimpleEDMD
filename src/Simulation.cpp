#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <cmath>
#include <algorithm>
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

void Simulation::addSphere(Vec3d pos, double radius){
    positions_.push_back(pos);
    radii_.push_back(radius);
    ++nSpheres_;
}

Vec3d Simulation::applyPeriodicBC(const Vec3d& vec)const{
    Vec3d retVec = vec;
    for(size_t i = 0; i < 3; ++i){
        retVec[i] -= boxSize_ * int(vec[i] * (2.0 / boxSize_));
    }
    return retVec;
}

//WARNING: Needs re-evaluation
bool Simulation::raySphereIntersection(double radius, const Vec3d& pos, const Vec3d& dir, double& t)const{
    double B   = dot(dir, pos);
    double det = B * B - dot(dir, dir) * (dot(pos, pos) - radius * radius);
    if(det < 0.0f) return false;

    double t0 = (B + sqrt(det)) / dot(dir, dir);
    double t1 = (B - sqrt(det)) / dot(dir, dir);
    bool retValue = false;
    if(t0 > 0.0){
        t = t0;
        //mIntersection = t * dir; //Intersection point
        retValue = true;
    }
    if(t1 < t && t1 > 0.0){
        t = t1;
        //mIntersection = t * dir;
        retValue = true;
    }
    return retValue;
}

CollisionEvent* Simulation::getCollisionEvent(size_t pA, size_t pB)const{
    Vec3d posA = positions_[pA],  posB = positions_[pB];
    Vec3d velA = velocities_[pA], velB = velocities_[pB];

    Vec3d dist   = applyPeriodicBC(posB - posA);
    Vec3d relVel = velB - velA;

    Time time(0.0);
    bool isCollision = raySphereIntersection(radii_[pA] + radii_[pB], dist, relVel, time);

    if(isCollision) return new CollisionEvent(time, pA, pB);
    else return nullptr;
}

//NOTE: For simplicity, for now we assume equal mass spheres
void Simulation::runCollisionEvent(const CollisionEvent& event){
    Time time = event.time_;
    size_t pA = event.pA;
    size_t pB = event.pB;

    positions_[pA] += velocities_[pA] * (time - times_[pA]);
    positions_[pB] += velocities_[pB] * (time - times_[pB]);

    times_[pA] = times_[pB] = time;

    auto temp = velocities_[pA];
    velocities_[pA] = velocities_[pB];
    velocities_[pB] = temp;

    auto assocsA = collisionGraph_->getAssociations(pA);
    collisionGraph_->clear(pA);
    for(auto eRef: assocsA) eventManager_.deleteEvent(eRef);
    auto assocsB = collisionGraph_->getAssociations(pB);
    collisionGraph_->clear(pB);
    for(auto eRef: assocsB) eventManager_.deleteEvent(eRef);

    //Recalculate collision events
    std::vector<CollisionEvent*> eventsA;
    std::vector<CollisionEvent*> eventsB;
    for(size_t n = 0; n < nSpheres_; ++n){
        if(n != pA && n != pB){
            auto eventA = getCollisionEvent(pA, n);
            auto eventB = getCollisionEvent(pB, n);
            if(eventA) eventsA.push_back(eventA);
            if(eventB) eventsB.push_back(eventB);
        }
        auto cmp = [](const CollisionEvent* a, const CollisionEvent* b){
            return (a->time_ < b->time_);
        };
        if(!eventsA.empty()){
            auto minEventA = *std::min_element(eventsA.begin(), eventsA.end(), cmp);
            EventRef refA = eventManager_.queueEvent(minEventA);
            collisionGraph_->addEdge(minEventA->pA, minEventA->pB, refA);
        }
        if(!eventsB.empty()){
            auto minEventB = *std::min_element(eventsB.begin(), eventsB.end(), cmp);
            EventRef refB = eventManager_.queueEvent(minEventB);
            collisionGraph_->addEdge(minEventB->pA, minEventB->pB, refB);
        }
    }
}

bool Simulation::init(void){
    if(nSpheres_) collisionGraph_ = new Graph(nSpheres_);
    else return false;

    //Initialize particle times to zero
    times_.resize(nSpheres_, 0.0);

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
    for(size_t i = 0; i < nSpheres_; ++i){
        std::vector<CollisionEvent*> events;
        for(size_t j = i + 1; j < nSpheres_; ++j){
            auto event = getCollisionEvent(i, j);
            if(event) events.push_back(event);
        }
        if(!events.empty()){
            CollisionEvent* earliestCollision = *std::min_element(events.begin(), events.end(),
                [](const CollisionEvent* a, const CollisionEvent* b){
                    return (a->time_ < b->time_);
                }
            );
            EventRef ref = eventManager_.queueEvent(earliestCollision);
            collisionGraph_->addEdge(earliestCollision->pA, earliestCollision->pB, ref);
        }
    }

    return true;
}

void Simulation::run(void){

    Time endTime = 100.0;
    
    bool running = true;
    while(running){
        const Event* nextEvent = eventManager_.getNextEvent();
        time_ = nextEvent->time_;
        std::cout << time_ << std::endl;

        if(time_ >= endTime) break;

        switch(nextEvent->getType()){
        case EVT_COLLISION:{
                auto collisionEvent = *static_cast<const CollisionEvent*>(nextEvent);
                runCollisionEvent(collisionEvent);
            }
            break;
        default:
            running = false;
            break;
        }
    }
}

