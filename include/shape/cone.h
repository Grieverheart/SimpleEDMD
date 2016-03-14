#ifndef EDMD_SHAPE_CONE_H
#define EDMD_SHAPE_CONE_H

#include "convex.h"

namespace shape{

    class Cone: public Convex{
    public:
        Cone(void){}

        explicit Cone(double base_radius, double height);

        friend void serialize(Archive& ar, const Cone&);
        friend void deserialize(Archive& ar, Cone*);

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
        double sintheta_;
        double in_radius_, out_radius_;
        double volume_;
    };

    inline clam::Vec3d Cone::support(const clam::Vec3d& dir)const{
        double test = dir[1] / dir.length();
        if(test >= sintheta_) return clam::Vec3d(0.0, half_height_, 0.0);
        else if(test < sintheta_ && (dir[0] != 0.0 || dir[2] != 0.0)){
            double length = sqrt(dir[0] * dir[0] + dir[2] * dir[2]);
            return clam::Vec3d(base_radius_ * dir[0] / length, -half_height_, base_radius_ * dir[2] / length);
        }
        else return clam::Vec3d(0.0, -half_height_, 0.0);
    }

}

#endif
