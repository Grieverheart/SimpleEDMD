#ifndef EDMD_BOUNDING_VOLUME_VARIANT_H
#define EDMD_BOUNDING_VOLUME_VARIANT_H

#include "shape/box.h"
#include "shape/sphere.h"
#include "bounding_volume_variant_fwd.h"
#include <boost/variant.hpp>

class BoundingVolumeOutRadiusVisitor: public boost::static_visitor<double>{
public:
    double operator()(const shape::Box& box)const{
        return box.out_radius();
    }

    double operator()(const shape::Sphere& sph)const{
        return sph.radius();
    }
};

inline double bv_outradius(const bounding_volume::Variant& shape){
    return boost::apply_visitor(BoundingVolumeOutRadiusVisitor(), shape);
}

#endif
