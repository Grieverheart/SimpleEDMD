#include "shape/cone.h"
#include "serialization/common.h"

namespace shape{
    //TODO: Move to separate header file
    template<typename T>
    static inline T sqr(T val){
        return val * val;
    }

    Cone::Cone(double base_radius, double height):
        base_radius_(base_radius), half_height_(0.5 * height),
        sintheta_(base_radius_ / sqrt(sqr(base_radius) + sqr(height))),
        in_radius_(sintheta_ * half_height_), out_radius_(sqrt(sqr(base_radius_) + sqr(half_height_))),
        volume_(M_PI * sqr(base_radius) * height / 3.0)
    {}

    void serialize(Archive& ar, const Cone& cone){
        serialize(ar, cone.base_radius_);
        serialize(ar, cone.half_height_);
        serialize(ar, cone.sintheta_);
        serialize(ar, cone.in_radius_);
        serialize(ar, cone.out_radius_);
        serialize(ar, cone.volume_);
    }

    void deserialize(Archive& ar, Cone* cone){
        deserialize(ar, &cone->base_radius_);
        deserialize(ar, &cone->half_height_);
        deserialize(ar, &cone->sintheta_);
        deserialize(ar, &cone->in_radius_);
        deserialize(ar, &cone->out_radius_);
        deserialize(ar, &cone->volume_);
    }

}

