#ifndef __PARTICLE_H
#define __PARTICLE_H

#include "clam.h"
#include "transformation.h"

struct Particle{
    Particle(void):
        time(0.0)
    {}
    double time;
    Transformation xform_;
    clam::Vec3d vel;
    clam::Vec3d ang_vel;
    size_t shape_id;
};

#endif
