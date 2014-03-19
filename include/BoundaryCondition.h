#ifndef __BOUNDARY_CONDITION_H
#define __BOUNDARY_CONDITION_H

#include "Vec.h"

class CubicPBC{
public:
    CubicPBC(void){}

    CubicPBC(double boxSize):
        boxSize_(boxSize), iBoxSize2_(2.0 / boxSize_)
    {}

    Vec3d minImage(const Vec3d& vec)const{
        Vec3d retVec;
    
        int k = vec.x * iBoxSize2_;
        retVec.x = vec.x - k * boxSize_;
        k = retVec.x * iBoxSize2_;
        retVec.x = retVec.x - k * boxSize_;
    
        int l = vec.y * iBoxSize2_;
        retVec.y = vec.y - l * boxSize_;
        l = retVec.y * iBoxSize2_;
        retVec.y = retVec.y - l * boxSize_;
    
        int m = vec.z * iBoxSize2_;
        retVec.z = vec.z - m * boxSize_;
        m = retVec.z * iBoxSize2_;
        retVec.z = retVec.z - m * boxSize_;
    
        return retVec;
    }

    Vec3d apply(const Vec3d& pos)const{
        Vec3d retVec;
        retVec.x = pos.x - int(pos.x * iBoxSize2_ - 1.0) * boxSize_;
        retVec.y = pos.y - int(pos.y * iBoxSize2_ - 1.0) * boxSize_;
        retVec.z = pos.z - int(pos.z * iBoxSize2_ - 1.0) * boxSize_;
        return retVec;
    }

    double getSize(void)const{
        return boxSize_;
    }

    void setSize(double boxSize){
        boxSize_   = boxSize;
        iBoxSize2_ = 2.0 / boxSize_;
    }

private:
    double boxSize_;
    double iBoxSize2_;
};

#endif
