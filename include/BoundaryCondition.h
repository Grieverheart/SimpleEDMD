#ifndef __BOUNDARY_CONDITION_H
#define __BOUNDARY_CONDITION_H

#include "clam.h"
#include "serialization/archive.h"
#include <cassert>
#include <cstdio>

//Assume y > 0.0
inline double wrap(double x, double y){
    if(x >= 0.0){
        if(x < y) return x;
        else if(x < 2.0 * y) return x - y;
    }
    else if(x >= -y) return x + y;

    // general case
    double m = x - y * int(x / y);
    // handle boundary cases resulting from floating-point limited accuracy:
    if(m >= y) return 0.0;
    if(m < 0.0){
        if(y == y + m) return 0.0;
        else return y + m;
    }

    return m;
}

class RectangularPBC{
public:
    RectangularPBC(void){}

    RectangularPBC(const clam::Vec3d& size):
        size_(size)
    {}

    clam::Vec3d minImage(const clam::Vec3d& vec)const{
        clam::Vec3d retVec;
        for(int i = 0; i < 3; ++i) retVec[i] = wrap(vec[i] + 0.5 * size_[i], size_[i]) - 0.5 * size_[i];
        return retVec;
    }

    clam::Vec3d apply(const clam::Vec3d& pos)const{
        clam::Vec3d retVec;
        for(int i = 0; i < 3; ++i) retVec[i] = wrap(pos[i], size_[i]);
        return retVec;
    }

    clam::Vec3d getSize(void)const{
        return size_;
    }

    void setSize(const clam::Vec3d& size){
        size_ = size;
    }

    void serialize(Archive& ar)const{
        size_.serialize(ar);
    }

private:
    clam::Vec3d size_;
};

#endif
