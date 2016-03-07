#include "overlap/obb.h"

bool overlap::obb_overlap(const Transform& box_a, const shape::Box& shape_a, const Transform& box_b, const shape::Box& shape_b, double margin){
    auto hsa = box_a.size_ * shape_a.extent() + margin;
    auto hsb = box_b.size_ * shape_b.extent() + margin;

    double e[9];
    box_b.rot_.to_matrix(e);

    double a[9] = {
        fabs(e[0]), fabs(e[3]), fabs(e[6]),
        fabs(e[1]), fabs(e[4]), fabs(e[7]),
        fabs(e[2]), fabs(e[5]), fabs(e[8])
    };

    const auto& pos = box_b.pos_;

    if(   (fabs(pos[0]) > hsa[0] + a[0] * hsb[0] + a[1] * hsb[1] + a[2] * hsb[2])
       || (fabs(pos[1]) > hsa[1] + a[3] * hsb[0] + a[4] * hsb[1] + a[5] * hsb[2])
       || (fabs(pos[2]) > hsa[2] + a[6] * hsb[0] + a[7] * hsb[1] + a[8] * hsb[2])
       || (fabs(e[0] * pos[0] + e[1] * pos[1] + e[2] * pos[2]) > hsb[0] + a[0] * hsa[0] + a[3] * hsa[1] + a[6] * hsa[2])
       || (fabs(e[3] * pos[0] + e[4] * pos[1] + e[5] * pos[2]) > hsb[1] + a[1] * hsa[0] + a[4] * hsa[1] + a[7] * hsa[2])
       || (fabs(e[6] * pos[0] + e[7] * pos[1] + e[8] * pos[2]) > hsb[2] + a[2] * hsa[0] + a[5] * hsa[1] + a[8] * hsa[2])
       || (fabs(pos[2] * e[1] - pos[1] * e[2]) > hsa[1] * a[6] + hsa[2] * a[3] + hsb[1] * a[2] + hsb[2] * a[1])
       || (fabs(pos[2] * e[4] - pos[1] * e[5]) > hsa[1] * a[7] + hsa[2] * a[4] + hsb[0] * a[2] + hsb[2] * a[0])
       || (fabs(pos[2] * e[7] - pos[1] * e[8]) > hsa[1] * a[8] + hsa[2] * a[5] + hsb[0] * a[1] + hsb[1] * a[0])
       || (fabs(pos[0] * e[2] - pos[2] * e[0]) > hsa[0] * a[6] + hsa[2] * a[0] + hsb[1] * a[5] + hsb[2] * a[4])
       || (fabs(pos[0] * e[5] - pos[2] * e[3]) > hsa[0] * a[7] + hsa[2] * a[1] + hsb[0] * a[5] + hsb[2] * a[3])
       || (fabs(pos[0] * e[8] - pos[2] * e[6]) > hsa[0] * a[8] + hsa[2] * a[2] + hsb[0] * a[4] + hsb[1] * a[3])
       || (fabs(pos[1] * e[0] - pos[0] * e[1]) > hsa[0] * a[3] + hsa[1] * a[0] + hsb[1] * a[8] + hsb[2] * a[7])
       || (fabs(pos[1] * e[3] - pos[0] * e[4]) > hsa[0] * a[4] + hsa[1] * a[1] + hsb[0] * a[8] + hsb[2] * a[6])
       || (fabs(pos[1] * e[6] - pos[0] * e[7]) > hsa[0] * a[5] + hsa[1] * a[2] + hsb[0] * a[7] + hsb[1] * a[6])
     ) return false;

    return true;
}
