#ifndef EDMD_SHAPE_VARIANT_FWD_H
#define EDMD_SHAPE_VARIANT_FWD_H

//@diary: For a generic shape, we can use a boost::variant. This solution is the most
//efficient and is similar to what we would have done in functional programming.
//The problem with the variant is that users cannot implement their own shapes without
//recompiling the library.

#include <boost/variant/variant_fwd.hpp>
#include <type_traits>

namespace shape{

    class Polyhedron;
    class Sphere;
    class Box;
    class Cone;
    class Cylinder;

    using Variant = boost::variant<Polyhedron, Sphere, Box, Cone, Cylinder>;
}

#endif
