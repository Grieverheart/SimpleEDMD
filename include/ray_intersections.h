#ifndef __RAY_INTERSECTIONS_H
#define __RAY_INTERSECTIONS_H

#include <cstdio>
#include <cmath>
#include "Vec.h"

inline bool raySphereIntersection(double radius, const Vec3d& pos, const Vec3d& dir, double& t){
    double dirInvLength = 1.0 / sqrt(dot(dir, dir));
    Vec3d dn  = dir * dirInvLength; //Normalize
    double s  = dot(pos, dn);
    double l2 = dot(pos, pos);
    double r2 = radius * radius;
    if(s < 0.0 && l2 > r2) return false;

    double m2 = l2 - s * s;
    if(m2 > r2) return false;

    double q = sqrt(r2 - m2);
    if(l2 > r2) t = s - q;
    else{
        t = s + q;
        printf("%f, %f\n", l2, r2);
    }
    t *= dirInvLength;

    return true;
}

//NOTE: Assumes that the ray is inside the box. Returns the new cell offset.
inline int rayCellIntersection(const Vec3d& cellSize, const Vec3d& rpos, const Vec3d& dir, double& t){
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
