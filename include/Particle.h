#ifndef __PARTICLE_H
#define __PARTICLE_H

#include "Vec.h"

struct Particle{
    Particle(void):
        time(0.0)
    {}
    double time;
    double radius;
    Vec3d  pos;
    Vec3d  vel;
};

#endif
