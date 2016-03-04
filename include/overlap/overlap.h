#ifndef EDMD_OVERLAP_OVERLAP_H
#define EDMD_OVERLAP_OVERLAP_H

#include "shape/variant_fwd.h"
#include "transform.h"

struct Particle;

namespace overlap{

    clam::Vec3d shape_distance(const Transform&, const shape::Variant&, const Transform&, const shape::Variant&);
    bool shape_overlap(const Transform&, const shape::Variant&, const Transform&, const shape::Variant&, double feather = 0.0);
    bool shape_raycast(const Transform&, const shape::Variant&, const Transform&, const shape::Variant&, const clam::Vec3d& ray_dir, double& distance, clam::Vec3d& normal);

}

#endif
