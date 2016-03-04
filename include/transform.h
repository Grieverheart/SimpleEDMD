#ifndef __TRANSFORM_H
#define __TRANSFORM_H

#include "clam.h"

struct Transform{
    clam::Vec3d pos_;
    clam::Quatd rot_;
    double size_;
};

#endif
