#ifndef __BOUNDARY_CONDITION_H
#define __BOUNDARY_CONDITION_H

#include "clam.h"
#include <cassert>
#include <cstdio>

class RectangularPBC{
public:
    RectangularPBC(void){}

    RectangularPBC(const clam::Vec3d& size):
        size_(size), isize_(1.0 / size_)
    {}

    //TODO: Check how accuracy affects simulation.

    //NOTE: Note sure if this always works.
    //Assumes -3.5B < Dx < Inf
    clam::Vec3d minImage(const clam::Vec3d& vec)const{
        clam::Vec3d retVec;
        for(int i = 0; i < 3; ++i){
            retVec[i] = vec[i] + (3 - int(vec[i] * isize_[i] + 3.5)) * size_[i];
            assert((retVec[i] - 0.5 * size_[i] < 1.0e-14) && (retVec[i] + 0.5 * size_[i] > -1.0e-14));
        }
    
        return retVec;
    }

    //Assumes -3.0B < pos < Inf
    clam::Vec3d apply(const clam::Vec3d& pos)const{
        clam::Vec3d retVec;
        for(int i = 0; i < 3; ++i){
            //NOTE: It is important to do the computations like this, for precision
            //reasons. If we were to elliminate terms, we would not get accurate results.
            //The same goes for using division instead of multiplying by isize_.
            double temp = pos[i] + 3.0 * size_[i];
            retVec[i] = temp - int(temp / size_[i]) * size_[i];
            assert((retVec[i] < size_[i]) && (retVec[i] >= 0.0));
        }

        return retVec;
    }

    clam::Vec3d getSize(void)const{
        return size_;
    }

    void setSize(const clam::Vec3d& size){
        size_  = size;
        isize_ = 1.0 / size_;
    }

private:
    clam::Vec3d size_;
    clam::Vec3d isize_;
};

#endif
