#ifndef EDMD_SHAPE_SPHERE_H
#define EDMD_SHAPE_SPHERE_H

#include "convex.h"
#include "serialization/archive.h"

namespace shape{

    class Sphere: public Convex{
    public:
        clam::Vec3d support(const clam::Vec3d& dir)const{
            return dir / dir.length();
        }

        double radius(void)const{
            return 1.0;
        }

        double volume(void)const{
            return (4.0 * M_PI / 3.0);
        }

    };

    inline void serialize(Archive& ar, const Sphere& sph){
    }

}

#endif
