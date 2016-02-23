#include "shape/box.h"
#include <algorithm>

namespace shape{

    Box::Box(const clam::Vec3d& dims):
        extent_(0.5 * dims),
        in_radius_(std::min(extent_[0], std::min(extent_[1], extent_[2]))), out_radius_(extent_.length()),
        volume_(dims[0] * dims[1] * dims[2])
    {}

}
