#ifndef __BOUNDARY_CONDITION_H
#define __BOUNDARY_CONDITION_H

#include "clam.h"

class CubicPBC{
public:
    CubicPBC(void){}

    CubicPBC(double boxSize):
        boxSize_(boxSize), iBoxSize2_(2.0 / boxSize_)
    {}

    clam::Vec3d minImage(const clam::Vec3d& vec)const{
        clam::Vec3d retVec;
    
        int k = vec[0] * iBoxSize2_;
        retVec[0] = vec[0] - k * boxSize_;
        k = retVec[0] * iBoxSize2_;
        retVec[0] = retVec[0] - k * boxSize_;
    
        int l = vec[1] * iBoxSize2_;
        retVec[1] = vec[1] - l * boxSize_;
        l = retVec[1] * iBoxSize2_;
        retVec[1] = retVec[1] - l * boxSize_;
    
        int m = vec[2] * iBoxSize2_;
        retVec[2] = vec[2] - m * boxSize_;
        m = retVec[2] * iBoxSize2_;
        retVec[2] = retVec[2] - m * boxSize_;
    
        return retVec;
    }

    clam::Vec3d apply(const clam::Vec3d& pos)const{
        clam::Vec3d retVec;
        //retVec[0] = pos[0] - int(pos[0] * iBoxSize2_ - 1.0) * boxSize_;
        //retVec[1] = pos[1] - int(pos[1] * iBoxSize2_ - 1.0) * boxSize_;
        //retVec[2] = pos[2] - int(pos[2] * iBoxSize2_ - 1.0) * boxSize_;
        for(int i = 0; i < 3; ++i){
            retVec[i] = pos[i];
            while(retVec[i] >= boxSize_) retVec[i] -= boxSize_;
            while(retVec[i] < 0.0) retVec[i] += boxSize_;
        }
        if(retVec[0] < 0.0 || retVec[1] < 0.0 || retVec[2] < 0.0){
            printf("%f, %f, %f\n", pos[0], pos[1], pos[2]);
        }
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
