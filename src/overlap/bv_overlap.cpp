#include "overlap/bv_overlap.h"
#include "bounding_volume_variant.h"
#include "overlap/gjk.h"

namespace overlap{

    namespace{
        //TODO: Move to separate header file
        template<typename T>
        static inline T sqr(T val){
            return val * val;
        }
    }

    class BoundingVolumeOverlapVisitor: public boost::static_visitor<bool> {
    public:
        BoundingVolumeOverlapVisitor(const Transform& pa, const Transform& pb, double feather):
            pa_(pa), pb_(pb), feather_(feather)
        {}

        //TODO: More efficient implementation. AABB vs Sphere.
        //Also probably not strictly correct as it checks for overlap between
        //rounded obb and sphere.
        template<typename T, typename U>
        bool operator()(const T& a, const U& b)const{
            return gjk_boolean(pa_, a, pb_, b, 2.0 * feather_);
        }

        bool operator()(const shape::Sphere& a, const shape::Sphere& b)const{
            return (pb_.pos_ - pa_.pos_).length2() < sqr(pa_.size_ * a.radius() + pb_.size_ * b.radius() + 2.0 * feather_);
        }

        bool operator()(const shape::Box& ba, const shape::Box& bb)const{
            auto hsa = pa_.size_ * ba.extent() + feather_;
            auto hsb = pb_.size_ * bb.extent() + feather_;

            double e[9];
            pb_.rot_.to_matrix(e);

            double a[9] = {
                fabs(e[0]), fabs(e[3]), fabs(e[6]),
                fabs(e[1]), fabs(e[4]), fabs(e[7]),
                fabs(e[2]), fabs(e[5]), fabs(e[8])
            };

            const auto& pos = pb_.pos_;

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

    private:
        const Transform& pa_;
        const Transform& pb_;
        double feather_;
    };

    bool bv_overlap(const Transform& ta, const bounding_volume::Variant& sa, const Transform& tb, const bounding_volume::Variant& sb, double feather){
        return boost::apply_visitor(BoundingVolumeOverlapVisitor(ta, tb, feather), sa, sb);
    }
}

