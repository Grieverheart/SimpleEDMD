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

            template<typename T>
            clam::Vec3d operator()(const shape::Complex& a, const T& b)const{
                return 0.0;
            }

            template<typename T>
            clam::Vec3d operator()(const T& a, const shape::Complex& b)const{
                return 0.0;
            }

            clam::Vec3d operator()(const shape::Complex& a, const shape::Complex& b)const{
                return 0.0;
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

            template<typename T>
            bool operator()(const shape::Complex& a, const T& b)const{
                return false;
            }

            template<typename T>
            bool operator()(const T& a, const shape::Complex& b)const{
                return false;
            }

            bool operator()(const shape::Complex& a, const shape::Complex& b)const{
                return false;
            }

        private:
            const Transformation& pa_;
            const Transformation& pb_;
            double feather_;
        };

        class ShapeRaycastVisitor: public boost::static_visitor<bool> {
        public:
            ShapeRaycastVisitor(const Transformation& pa, const Transformation& pb, const clam::Vec3d& ray_dir, double& distance, clam::Vec3d& normal):
                pa_(pa), pb_(pb), ray_dir_(ray_dir), dist_(distance), normal_(normal)
            {}

            template<typename T, typename U>
            bool operator()(const T& a, const U& b)const{
                return gjk_raycast(pa_, a, pb_, b, ray_dir_, dist_, normal_);
            }

            bool operator()(const shape::Sphere& a, const shape::Sphere& b)const{
                return sphere_raycast(pa_.size_ * a.radius() + pb_.size_ * b.radius(), pa_.pos_ - pb_.pos_, ray_dir_, dist_, &normal_);
            }

            template<typename T>
            bool operator()(const shape::Complex& a, const T& b)const{
                return false;
            }

            template<typename T>
            bool operator()(const T& a, const shape::Complex& b)const{
                return false;
            }

            bool operator()(const shape::Complex& a, const shape::Complex& b)const{
                return false;
            }

        private:
            const Transformation& pa_;
            const Transformation& pb_;
            const clam::Vec3d& ray_dir_;
            double& dist_;
            clam::Vec3d& normal_;
        };
    }//namespace detail

    clam::Vec3d shape_distance(const Transformation& pa, const shape::Variant& a, const Transformation& pb, const shape::Variant& b){
        return boost::apply_visitor(ShapeDistanceVisitor(pa, pb), a, b);
    }

    bool shape_overlap(const Transformation& pa, const shape::Variant& a, const Transformation& pb, const shape::Variant& b, double feather){
        return boost::apply_visitor(ShapeOverlapVisitor(pa, pb, feather), a, b);
    }

    bool shape_raycast(const Transformation& pa, const shape::Variant& a, const Transformation& pb, const shape::Variant& b, const clam::Vec3d& ray_dir, double& distance, clam::Vec3d& normal){
        return boost::apply_visitor(ShapeRaycastVisitor(pa, pb, ray_dir, distance, normal), a, b);
    }

}
