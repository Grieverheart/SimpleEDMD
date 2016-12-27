#ifndef EDMD_SHAPE_LEAF_CYLINDER_H
#define EDMD_SHAPE_LEAF_CYLINDER_H

#include "convex.h"
#include <cmath>

namespace shape{

    //NOTE: w < l!
    /* Vesica Piscis */
    class LeafCylinder: public Convex{
    public:
        LeafCylinder(void){}

        explicit LeafCylinder(double width, double length, double height);

        friend void serialize(Archive& ar, const LeafCylinder&);
        friend void deserialize(Archive& ar, LeafCylinder*);

        double in_radius(void)const{
            return in_radius_;
        }
        double out_radius(void)const{
            return out_radius_;
        }
        double volume(void)const{
            return volume_;
        }
        double width(void)const{
            return 2.0 * half_width_;
        }
        double length(void)const{
            return 2.0 * half_length_;
        }
        double height(void)const{
            return 2.0 * half_height_;
        }

        clam::Vec3d support(const clam::Vec3d&)const;
    private:
        double half_width_, half_length_, half_height_;
        double circle_radius_, circle_distance_;
        double in_radius_, out_radius_;
        double volume_;
    };

    inline clam::Vec3d LeafCylinder::support(const clam::Vec3d& dir)const{
        double x = 0.0, z = 0.0;
        if(dir[0] != 0.0 || dir[2] != 0.0){
            double l = sqrt(dir[0] * dir[0] + dir[2] * dir[2]);
            double test = dir[2] / l;
            if(test >= half_length_ / circle_radius_) z = half_length_;
            else if(test <= -half_length_ / circle_radius_) z = -half_length_;
            else{
                x = circle_radius_ * dir[0] / l - std::copysign(circle_distance_, dir[0]);
                z = circle_radius_ * dir[2] / l;
            }
        }

        return clam::Vec3d(x, std::copysign(half_height_, dir[1]), z);
    }

}

#endif
