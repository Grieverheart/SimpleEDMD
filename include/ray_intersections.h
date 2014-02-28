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

#endif
