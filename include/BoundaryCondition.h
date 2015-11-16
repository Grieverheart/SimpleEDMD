#ifndef __BOUNDARY_CONDITION_H
#define __BOUNDARY_CONDITION_H

#include "clam.h"
#include <cassert>

class RectangularPBC{
public:
    RectangularPBC(void){}

    RectangularPBC(const clam::Vec3d& size):
        size_(size), isize_(2.0 / size_)
    {}

    //Assum -1.5B < Dx < 1.5B
    clam::Vec3d minImage(const clam::Vec3d& vec)const{
        clam::Vec3d retVec;
    
        retVec[0] = vec[0] - int(vec[0] * isize_[0]) * size_[0];
        retVec[1] = vec[1] - int(vec[1] * isize_[1]) * size_[1];
        retVec[2] = vec[2] - int(vec[2] * isize_[2]) * size_[2];

        retVec[0] -= int(retVec[0] * isize_[0]) * size_[0];
        retVec[1] -= int(retVec[1] * isize_[1]) * size_[1];
        retVec[2] -= int(retVec[2] * isize_[2]) * size_[2];

        assert((retVec[0] < 0.5 * size_[0]) && (retVec[0] > -0.5 * size_[0]));
        assert((retVec[1] < 0.5 * size_[1]) && (retVec[1] > -0.5 * size_[1]));
        assert((retVec[2] < 0.5 * size_[2]) && (retVec[2] > -0.5 * size_[2]));
    
        return retVec;
    }

    //TODO: Make this faster.
    clam::Vec3d apply(const clam::Vec3d& pos)const{
        clam::Vec3d retVec;

        for(int i = 0; i < 3; ++i){
            retVec[i] = pos[i];
            while(retVec[i] < 0.0) retVec[i] += size_[i];
            while(retVec[i] >= size_[i]) retVec[i] -= size_[i];
        }

        assert((retVec[0] < size_[0]) && (retVec[0] >= 0.0));
        assert((retVec[1] < size_[1]) && (retVec[1] >= 0.0));
        assert((retVec[2] < size_[2]) && (retVec[2] >= 0.0));

        return retVec;
    }

    clam::Vec3d getSize(void)const{
        return size_;
    }

    void setSize(const clam::Vec3d& size){
        size_  = size;
        isize_ = 2.0 / size_;
    }

private:
    clam::Vec3d size_;
    clam::Vec3d isize_;
};

#endif
