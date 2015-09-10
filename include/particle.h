#ifndef __PARTICLE_H
#define __PARTICLE_H

#include "clam.h"

struct Particle{
    Particle(void):
        time(0.0)
    {}
    double time;
    double size;
    clam::Vec3d pos;
    clam::Quatd rot;
    clam::Vec3d vel;
    size_t shape_id;
};

#endif
