#ifndef __TRANSFORMATION_H
#define __TRANSFORMATION_H

#include "clam.h"

//TODO: Add a multiplication operator for combining transformations.
struct Transformation{
    Transformation(void){}

    Transformation(const clam::Vec3d& pos, const clam::Quatd& rot, double size = 1.0):
        pos_(pos), rot_(rot), size_(size)
    {}

    clam::Vec3d pos_;
    clam::Quatd rot_;
    double size_;
};

#endif
