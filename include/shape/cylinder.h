#ifndef EDMD_SHAPE_CYLINDER_H
#define EDMD_SHAPE_CYLINDER_H

#include "convex.h"
#include <cmath>

namespace shape{

    class Cylinder: public Convex{
    public:
        Cylinder(void){}

        explicit Cylinder(double base_radius, double height);

        friend void serialize(Archive& ar, const Cylinder&);
        friend void deserialize(Archive& ar, Cylinder*);

        double in_radius(void)const{
            return in_radius_;
        }
        double out_radius(void)const{
            return out_radius_;
        }
        double volume(void)const{
            return volume_;
        }
        double base_radius(void)const{
            return base_radius_;
        }
        double height(void)const{
            return 2.0 * half_height_;
        }

        clam::Vec3d support(const clam::Vec3d&)const;
    private:
        double base_radius_, half_height_;
        double in_radius_, out_radius_;
        double volume_;
    };

    inline clam::Vec3d Cylinder::support(const clam::Vec3d& dir)const{
        if(dir[0] != 0.0 || dir[2] != 0.0){
            double length = sqrt(dir[0] * dir[0] + dir[2] * dir[2]);
            return clam::Vec3d(base_radius_ * dir[0] / length, std::copysign(half_height_, dir[1]), base_radius_ * dir[2] / length);
        }
        else return clam::Vec3d(0.0, std::copysign(half_height_, dir[1]), 0.0);
    }

}

#endif
