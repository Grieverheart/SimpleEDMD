#include "overlap/overlap.h"
#include "particle.h"
#include "shape/variant.h"
#include "overlap/gjk.h"
#include "overlap/ray_casting.h"

namespace overlap{

    namespace{

        class ShapeDistanceVisitor: public boost::static_visitor<clam::Vec3d> {
        public:
            ShapeDistanceVisitor(const Transformation& pa, const Transformation& pb):
                pa_(pa), pb_(pb)
            {}

            template<typename T, typename U>
            clam::Vec3d operator()(const T& a, const U& b)const{
                return gjk_distance(pa_, a, pb_, b);
            }

            clam::Vec3d operator()(const shape::Sphere& a, const shape::Sphere& b)const{
                return pb_.pos_ - pa_.pos_ - pa_.size_ * a.radius() - pb_.size_ * b.radius();
            }

            clam::Vec3d operator()(const shape::Complex& a, const shape::Complex& b)const{
                clam::Vec3d min_dist = clam::Vec3d(std::numeric_limits<double>::max());
                for(auto subshape_a: a.shapes()){
                    auto xform_a = pa_ * subshape_a.xform_;
                    for(auto subshape_b: b.shapes()){
                        auto xform_b = pb_ * subshape_b.xform_;
                        clam::Vec3d dist = shape_distance(xform_a, *subshape_a.shape_, xform_b, *subshape_b.shape_);
                        if(dist.length2() < min_dist.length2()) min_dist = dist;
                    }
                }
                return min_dist;
            }

            template<typename T>
            clam::Vec3d operator()(const shape::Complex& a, const T& b)const{
                clam::Vec3d min_dist = clam::Vec3d(std::numeric_limits<double>::max());
                auto b_variant = shape::Variant(b);
                for(auto subshape_a: a.shapes()){
                    auto xform_a = pa_ * subshape_a.xform_;
                    clam::Vec3d dist = shape_distance(xform_a, *subshape_a.shape_, pb_, b_variant);
                    if(dist.length2() < min_dist.length2()) min_dist = dist;
                }
                return min_dist;
            }

            template<typename T>
            clam::Vec3d operator()(const T& a, const shape::Complex& b)const{
                clam::Vec3d min_dist = clam::Vec3d(std::numeric_limits<double>::max());
                auto a_variant = shape::Variant(a);
                for(auto subshape_b: b.shapes()){
                    auto xform_b = pb_ * subshape_b.xform_;
                    clam::Vec3d dist = shape_distance(pa_, a_variant, xform_b, *subshape_b.shape_);
                    if(dist.length2() < min_dist.length2()) min_dist = dist;
                }
                return min_dist;
            }

        private:
            const Transformation& pa_;
            const Transformation& pb_;
        };

        //TODO: Move to separate header file
        template<typename T>
        static inline T sqr(T val){
            return val * val;
        }

        class ShapeOverlapVisitor: public boost::static_visitor<bool> {
        public:
            ShapeOverlapVisitor(const Transformation& pa, const Transformation& pb, double feather):
                pa_(pa), pb_(pb), feather_(feather)
            {}

            template<typename T, typename U>
            bool operator()(const T& a, const U& b)const{
                return gjk_boolean(pa_, a, pb_, b, feather_);
            }

            bool operator()(const shape::Sphere& a, const shape::Sphere& b)const{
                return (pb_.pos_ - pa_.pos_).length2() < sqr(pa_.size_ * a.radius() + pb_.size_ * b.radius() + feather_);
            }

            bool operator()(const shape::Complex& a, const shape::Complex& b)const{
                for(auto subshape_a: a.shapes()){
                    auto xform_a = pa_ * subshape_a.xform_;
                    for(auto subshape_b: b.shapes()){
                        auto xform_b = pb_ * subshape_b.xform_;
                        if(shape_overlap(xform_a, *subshape_a.shape_, xform_b, *subshape_b.shape_, feather_)) return true;
                    }
                }
                return false;
            }

            template<typename T>
            bool operator()(const shape::Complex& a, const T& b)const{
                auto b_variant = shape::Variant(b);
                for(auto subshape_a: a.shapes()){
                    auto xform_a = pa_ * subshape_a.xform_;
                    if(shape_overlap(xform_a, *subshape_a.shape_, pb_, b_variant, feather_)) return true;
                }
                return false;
            }

            template<typename T>
            bool operator()(const T& a, const shape::Complex& b)const{
                auto a_variant = shape::Variant(a);
                for(auto subshape_b: b.shapes()){
                    auto xform_b = pb_ * subshape_b.xform_;
                    if(shape_overlap(pa_, a_variant, xform_b, *subshape_b.shape_, feather_)) return true;
                }
                return false;
            }

        private:
            const Transformation& pa_;
            const Transformation& pb_;
            double feather_;
        };
    }//namespace detail

    clam::Vec3d shape_distance(const Transformation& pa, const shape::Variant& a, const Transformation& pb, const shape::Variant& b){
        return boost::apply_visitor(ShapeDistanceVisitor(pa, pb), a, b);
    }

    bool shape_overlap(const Transformation& pa, const shape::Variant& a, const Transformation& pb, const shape::Variant& b, double feather){
        return boost::apply_visitor(ShapeOverlapVisitor(pa, pb, feather), a, b);
    }

}
