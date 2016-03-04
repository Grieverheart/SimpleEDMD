#ifndef __PARTICLE_H
#define __PARTICLE_H

#include "transform.h"

struct Particle{
    Particle(void):
        time(0.0)
    {}

    double time;
    Transform xform;
    clam::Vec3d vel;
    clam::Vec3d ang_vel;
    size_t shape_id;
};

#endif
