#include "shape/leaf_cylinder.h"
#include "serialization/common.h"
#include <algorithm>

namespace shape{
    //TODO: Move to separate header file
    template<typename T>
    static inline T sqr(T val){
        return val * val;
    }

    LeafCylinder::LeafCylinder(double width, double length, double height):
        half_width_(0.5 * width), half_length_(0.5 * length), half_height_(0.5 * height),
        circle_radius_(0.25 * (sqr(length) + sqr(width)) / width), circle_distance_(0.25 * (sqr(length) - sqr(width)) / width),
        in_radius_(std::min(half_height_, half_width_)), out_radius_(sqrt(sqr(half_height_) + sqr(half_length_)))
    {
        double theta = 2.0 * asin(half_length_ / circle_radius_);
        volume_ = height * sqr(circle_radius_) * (theta - sin(theta));
    }

    void serialize(Archive& ar, const LeafCylinder& leaf){
        serialize(ar, leaf.half_width_);
        serialize(ar, leaf.half_length_);
        serialize(ar, leaf.half_height_);
        serialize(ar, leaf.circle_radius_);
        serialize(ar, leaf.circle_distance_);
        serialize(ar, leaf.in_radius_);
        serialize(ar, leaf.out_radius_);
        serialize(ar, leaf.volume_);
    }

    void deserialize(Archive& ar, LeafCylinder* leaf){
        deserialize(ar, &leaf->half_width_);
        deserialize(ar, &leaf->half_length_);
        deserialize(ar, &leaf->half_height_);
        deserialize(ar, &leaf->circle_radius_);
        deserialize(ar, &leaf->circle_distance_);
        deserialize(ar, &leaf->in_radius_);
        deserialize(ar, &leaf->out_radius_);
        deserialize(ar, &leaf->volume_);
    }


}
