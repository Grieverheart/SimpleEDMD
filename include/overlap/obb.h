#ifndef EDMD_OVERLAP_OBB_H
#define EDMD_OVERLAP_OBB_H

#include "shape/box.h"
#include "transform.h"

namespace overlap{
    //OBB-OBB overlap check. box_b's transformation should be with respect to box_a.
    //i.e. box_a is considered to be at the origin and not rotated.
    bool obb_overlap(const Transform& box_a, const shape::Box& shape_a, const Transform& box_b, const shape::Box& shape_b, double margin = 0.0);
}

#endif
