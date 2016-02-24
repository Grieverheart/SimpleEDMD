#ifndef EDMD_OVERLAP_OBB_H
#define EDMD_OVERLAP_OBB_H

#include "shape/box.h"
#include "bounding_box.h"

namespace overlap{
    bool obb_overlap(const BoundingBox& box_a, const shape::Box& shape_a, const BoundingBox& box_b, const shape::Box& shape_b, double margin = 0.0);
}

#endif
