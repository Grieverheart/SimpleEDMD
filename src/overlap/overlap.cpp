#include "overlap/overlap.h"
#include "particle.h"
#include "shape/variant.h"
#include "overlap/gjk.h"
#include "overlap/ray_casting.h"

namespace overlap{

    namespace{

        class ShapeDistanceVisitor: public boost::static_visitor<double> {
        public:
            ShapeDistanceVisitor(const Particle& pa, const Particle& pb):
                pa_(pa), pb_(pb)
            {}

            template<typename T, typename U>
            double operator()(const T& a, const U& b)const{
                return gjk_distance(pa_, a, pb_, b);
            }

        private:
            const Particle& pa_;
            const Particle& pb_;
        };

        template<>
        inline double ShapeDistanceVisitor::operator()(const shape::Sphere& a, const shape::Sphere& b)const{
            return (pb_.pos - pa_.pos).length() - pa_.size * a.radius() - pb_.size * b.radius();
        }

        //TODO: Move to separate header file
        template<typename T>
        static inline T sqr(T val){
            return val * val;
        }

        class ShapeOverlapVisitor: public boost::static_visitor<bool> {
        public:
            ShapeOverlapVisitor(const Particle& pa, const Particle& pb, double feather):
                pa_(pa), pb_(pb), feather_(feather)
            {}

            template<typename T, typename U>
            bool operator()(const T& a, const U& b)const{
                return gjk_boolean(pa_, a, pb_, b, feather_);
            }

        private:
            const Particle& pa_;
            const Particle& pb_;
            double feather_;
        };

        template<>
        inline bool ShapeOverlapVisitor::operator()(const shape::Sphere& a, const shape::Sphere& b)const{
            return (pb_.pos - pa_.pos).length2() < sqr(pa_.size * a.radius() + pb_.size * b.radius() + feather_);
        }

        class ShapeRaycastVisitor: public boost::static_visitor<bool> {
        public:
            ShapeRaycastVisitor(const Particle& pa, const Particle& pb, const clam::Vec3d& ray_dir, double& distance):
                pa_(pa), pb_(pb), ray_dir_(ray_dir), dist_(distance)
            {}

            template<typename T, typename U>
            bool operator()(const T& a, const U& b)const{
                return gjk_raycast(pa_, a, pb_, b, ray_dir_, dist_);
            }

        private:
            const Particle& pa_;
            const Particle& pb_;
            const clam::Vec3d& ray_dir_;
            double& dist_;
        };

        template<>
        inline bool ShapeRaycastVisitor::operator()(const shape::Sphere& a, const shape::Sphere& b)const{
            return sphere_raycast(pa_.size * a.radius() + pb_.size * b.radius(), pa_.pos - pb_.pos, ray_dir_, dist_);
        }
    }//namespace detail

    double shape_distance(const Particle& pa, const shape::Variant& a, const Particle& pb, const shape::Variant& b){
        return boost::apply_visitor(ShapeDistanceVisitor(pa, pb), a, b);
    }

    bool shape_overlap(const Particle& pa, const shape::Variant& a, const Particle& pb, const shape::Variant& b, double feather){
        return boost::apply_visitor(ShapeOverlapVisitor(pa, pb, feather), a, b);
    }

    bool shape_raycast(const Particle& pa, const shape::Variant& a, const Particle& pb, const shape::Variant& b, const clam::Vec3d& ray_dir, double& distance){
        return boost::apply_visitor(ShapeRaycastVisitor(pa, pb, ray_dir, distance), a, b);
    }

}
