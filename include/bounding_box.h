#ifndef __BOUNDING_BOX_H
#define __BOUNDING_BOX_H

#include "clam.h"

//TODO: Rename to Coordinates, and add 'size'
struct BoundingBox{
    clam::Vec3d pos_;
    clam::Quatd rot_;
};

#endif
