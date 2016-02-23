#ifndef MCH_SHAPE_BOX_H
#define MCH_SHAPE_BOX_H

#include "convex.h"
#include <cmath>

namespace shape{

    class Box: public Convex{
    public:
        Box(const clam::Vec3d& dims);

        double in_radius(void)const{
            return in_radius_;
        }
        double out_radius(void)const{
            return out_radius_;
        }
        double volume(void)const{
            return volume_;
        }
        clam::Vec3d extent(void)const{
            return extent_;
        }

        clam::Vec3d support(const clam::Vec3d&)const;
    private:
        clam::Vec3d extent_;
        double in_radius_, out_radius_;
        double volume_;
    };

    inline clam::Vec3d Box::support(const clam::Vec3d& dir)const{
        return clam::Vec3d(
            std::copysign(extent_[0], dir[0]),
            std::copysign(extent_[1], dir[1]),
            std::copysign(extent_[2], dir[2])
        );
    }

}

#endif
