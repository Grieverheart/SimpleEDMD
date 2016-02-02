#ifndef __PARTICLE_H
#define __PARTICLE_H

#include "serialization/archive.h"
#include "clam.h"

struct Particle{
    Particle(void):
        time(0.0)
    {}

    void serialize(Archive& ar){
        ar.write(&time, sizeof(double));
        ar.write(&size, sizeof(double));
        pos.serialize(ar);
        rot.serialize(ar);
        vel.serialize(ar);
        ang_vel.serialize(ar);
        ar.write(&shape_id, sizeof(double));
    }

    double time;
    double size;
    clam::Vec3d pos;
    clam::Quatd rot;
    clam::Vec3d vel;
    clam::Vec3d ang_vel;
    size_t shape_id;
};

#endif
