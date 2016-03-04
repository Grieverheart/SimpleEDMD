#ifndef EDMD_BOUNDING_VOLUME_VARIANT_FWD_H
#define EDMD_BOUNDING_VOLUME_VARIANT_FWD_H

#include <boost/variant/variant_fwd.hpp>

namespace shape{
    class Box;
    class Sphere;
}

namespace bounding_volume{

    enum eShapeTypes{
        BOX    = 0,
        SPHERE = 1,
        N_SHAPE_TYPES
    };

    using Variant = boost::variant<shape::Box, shape::Sphere>;
}

#endif
