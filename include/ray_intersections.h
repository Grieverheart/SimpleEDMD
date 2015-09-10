#ifndef __RAY_INTERSECTIONS_H
#define __RAY_INTERSECTIONS_H

#include <cstdio>
#include <cmath>
#include "clam.h"

inline bool raySphereIntersection(double radius, const clam::Vec3d& pos, const clam::Vec3d& dir, double& t){
    double s  = clam::dot(pos, dir);
    double l2 = clam::dot(pos, pos);
    double r2 = radius * radius;
    //NOTE: maybe we should just check s < 0.0. This means they are moving appart.
    //We only want to know if they are approaching eachother and l2 < r2.
    if(s < 0.0 && l2 > r2) return false;

    double idnorm = 1.0 / sqrt(clam::dot(dir, dir));
    s *= idnorm;
    double m2 = l2 - s * s;
    if(m2 > r2) return false;

    double q = sqrt(r2 - m2);
    if(__builtin_expect(l2 > r2, 1)) t = s - q;
    else{
        t = s + q;
        //printf("%f, %f\n", l2, r2);
    }
    t *= idnorm;

    return true;
}

//NOTE: Expects spheres to not be penetrating
inline bool raySphereIntersectionF(double radius, const clam::Vec3d& pos, const clam::Vec3d& dir, double& t){
    double s  = clam::dot(pos, dir);
    if(s < 0.0) return false;

    double l2 = clam::dot(pos, pos);
    double r2 = radius * radius;
    double idnorm = 1.0 / sqrt(clam::dot(dir, dir));
    s *= idnorm;
    double m2 = l2 - s * s;
    if(m2 > r2) return false;

    t = (s - sqrt(r2 - m2)) * idnorm;
    return true;
}

//NOTE: Assumes that the ray is inside the box. Returns the new cell offset.
inline int rayCellIntersection(const clam::Vec3d& cellSize, const clam::Vec3d& rpos, const clam::Vec3d& dir, double& t){
    //NOTE: Check for +-0
    t = (dir[0] < 0.0)? -rpos[0] / dir[0]: (cellSize[0] - rpos[0]) / dir[0];
    int cellOffset = !(dir[0] < 0.0);
    for(int i = 1; i < 3; ++i){
        bool isNegative = (dir[i] < 0.0);
        double dt = isNegative? -rpos[i] / dir[i]: (cellSize[i] - rpos[i]) / dir[i];
        if(dt < t){
            t = dt;
            cellOffset = 2 * i + !isNegative;
        }
    }
    return cellOffset;
}

#endif
