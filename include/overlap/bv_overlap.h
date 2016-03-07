#ifndef EDMD_OVERLAP_BV_OVERLAP_H
#define EDMD_OVERLAP_BV_OVERLAP_H

#include "bounding_volume_variant_fwd.h"
#include "transform.h"

struct Particle;

namespace overlap{

    bool bv_overlap(const Transform&, const bounding_volume::Variant&, const Transform&, const bounding_volume::Variant&, double feather = 0.0);

}

#endif
