#include <cstdio>
#include <cstring>
#include <cstdlib>
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

void Simulation::run(void){
    EventRef impendingCollisions[nSpheres_];
    EventRef impendingTransfers[nSpheres_];

    /* Pseudocode for 'run' function */
    //Initialize paricle velocities
    
    for(int i = 0; i < nSpheres_; ++i){
        for(int j = i + 1; j < nSpheres_; ++j){
            //get collision info
            //if collision, then store collision info
        }
        //Calculate when particle will move out of box
        //queue earliest collision event and transfer event
    }

    bool running = true;
    while(running){
        const Event* nextEvent = eventManager_.getNextEvent();
        switch(nextEvent->getType()){
        case EVT_COLLISION:{
                Time time = nextEvent->time_;
                auto collisionEvent = *static_cast<const CollisionEvent*>(nextEvent);
                times_[collisionEvent.pA] = times_[collisionEvent.pB] = time;
                //Change participating particles' velocity, time and position
                //Recalculate collision events
            }
            break;
        case EVT_TRANSFER:
            //Apply periodic boundary conditions for the particle and recalculate
            //collision events
            break;
        case EVT_OUTPUT:
            //progress all particles to current time and output configuration
            //schedule next output event
            break;
        case EVT_ENDSIM:
            running = false;
            break;
        }
    }
}

