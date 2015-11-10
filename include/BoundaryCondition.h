#ifndef __BOUNDARY_CONDITION_H
#define __BOUNDARY_CONDITION_H

#include "clam.h"
#include <cassert>

class CubicPBC{
public:
    CubicPBC(void){}

    CubicPBC(double boxSize):
        boxSize_(boxSize), iBoxSize2_(2.0 / boxSize_)
    {}

    //Assum -1.5B < Dx < 1.5B
    clam::Vec3d minImage(const clam::Vec3d& vec)const{
        clam::Vec3d retVec;
    
        retVec[0] = vec[0] - int(vec[0] * iBoxSize2_) * boxSize_;
        retVec[1] = vec[1] - int(vec[1] * iBoxSize2_) * boxSize_;
        retVec[2] = vec[2] - int(vec[2] * iBoxSize2_) * boxSize_;

        retVec[0] -= int(retVec[0] * iBoxSize2_) * boxSize_;
        retVec[1] -= int(retVec[1] * iBoxSize2_) * boxSize_;
        retVec[2] -= int(retVec[2] * iBoxSize2_) * boxSize_;

        assert((retVec[0] < 0.5 * boxSize_) && (retVec[0] > -0.5 * boxSize_));
        assert((retVec[1] < 0.5 * boxSize_) && (retVec[1] > -0.5 * boxSize_));
        assert((retVec[2] < 0.5 * boxSize_) && (retVec[2] > -0.5 * boxSize_));
    
        return retVec;
    }

    //TODO: Make this faster.
    clam::Vec3d apply(const clam::Vec3d& pos)const{
        clam::Vec3d retVec;
        //retVec[0] = pos[0] - int(pos[0] * iBoxSize2_ - 1.0) * boxSize_;
        //retVec[1] = pos[1] - int(pos[1] * iBoxSize2_ - 1.0) * boxSize_;
        //retVec[2] = pos[2] - int(pos[2] * iBoxSize2_ - 1.0) * boxSize_;
        for(int i = 0; i < 3; ++i){
            retVec[i] = pos[i];
            while(retVec[i] < 0.0) retVec[i] += boxSize_;
            while(retVec[i] >= boxSize_) retVec[i] -= boxSize_;
        }

        assert((retVec[0] < boxSize_) && (retVec[0] >= 0.0));
        assert((retVec[1] < boxSize_) && (retVec[1] >= 0.0));
        assert((retVec[2] < boxSize_) && (retVec[2] >= 0.0));

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
