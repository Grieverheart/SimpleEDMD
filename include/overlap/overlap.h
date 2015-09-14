#ifndef EDMD_OVERLAP_OVERLAP_H
#define EDMD_OVERLAP_OVERLAP_H

#include "shape/variant_fwd.h"
#include "clam.h"

struct Particle;

namespace overlap{

    double shape_distance(const Particle&, const shape::Variant&, const Particle&, const shape::Variant&);
    bool shape_overlap(const Particle&, const shape::Variant&, const Particle&, const shape::Variant&, double feather = 0.0);
    bool shape_raycast(const Particle&, const shape::Variant&, const Particle&, const shape::Variant&, const clam::Vec3d& ray_dir, double& distance, clam::Vec3d& normal);

}

#endif
