#ifndef EDMD_OVERLAP_OVERLAP_H
#define EDMD_OVERLAP_OVERLAP_H

#include "shape/variant_fwd.h"
#include "transformation.h"
#include "clam.h"

namespace overlap{

    //TODO: Perhaps it's better if we pass pos, rot, instead of the whole particle.
    clam::Vec3d shape_distance(const Transformation&, const shape::Variant&, const Transformation&, const shape::Variant&);
    bool shape_overlap(const Transformation&, const shape::Variant&, const Transformation&, const shape::Variant&, double feather = 0.0);

}

#endif
