#include "shape/complex.h"
#include "shape/variant.h"

namespace shape{

void Complex::update_outradius(void){
    const auto& shape_def = shapes_.back();
    double max_dist = shape_def.xform_.pos_.length() + shape_def.xform_.size_ * shape_outradius(*shape_def.shape_);
    if(max_dist > out_radius_) out_radius_ = max_dist;
}

}
