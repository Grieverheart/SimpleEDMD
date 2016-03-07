#include "shape/box.h"
#include "serialization/common.h"
#include <algorithm>

namespace shape{

    Box::Box(const clam::Vec3d& dims):
        extent_(0.5 * dims),
        in_radius_(std::min(extent_[0], std::min(extent_[1], extent_[2]))), out_radius_(extent_.length()),
        volume_(dims[0] * dims[1] * dims[2])
    {}

    void serialize(Archive& ar, const Box& box){
        serialize(ar, box.extent_);
        serialize(ar, box.in_radius_);
        serialize(ar, box.out_radius_);
        serialize(ar, box.volume_);
    }

    void deserialize(Archive& ar, Box* box){
        serialize(ar, &box->extent_);
        serialize(ar, &box->in_radius_);
        serialize(ar, &box->out_radius_);
        serialize(ar, &box->volume_);
    }

}
