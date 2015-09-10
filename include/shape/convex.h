#ifndef EDMD_SHAPE_CONVEX_H
#define EDMD_SHAPE_CONVEX_H

#include "clam.h"

namespace shape{

    class Convex{
    public:
        virtual ~Convex(void){};
        virtual clam::Vec3d support(const clam::Vec3d&)const = 0;
        virtual double volume(void)const = 0;
    };

}

#endif
