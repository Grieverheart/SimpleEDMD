#ifndef EDMD_OVERLAP_OBB_H
#define EDMD_OVERLAP_OBB_H

struct BBShape;
struct BoundingBox;

namespace overlap{
    bool obb_overlap(const BoundingBox& box_a, const BBShape& shape_a, const BoundingBox& box_b, const BBShape& shape_b);
}

#endif
