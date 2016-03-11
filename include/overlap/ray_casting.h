#ifndef EDMD_OVERLAP_RAY_CASTING_H
#define EDMD_OVERLAP_RAY_CASTING_H

namespace overlap{
    //If spheres are penetrating, returns false.
    inline bool sphere_raycast(double radius, const clam::Vec3d& pos, const clam::Vec3d& dir, double& t, clam::Vec3d* normal = nullptr){
        double s  = dot(pos, dir);
        double l2 = dot(pos, pos);
        double r2 = radius * radius;

        if(s < 0.0/* && l2 > r2*/) return false;

        double idnorm = 1.0 / dir.length();
        s *= idnorm;
        double m2 = l2 - s * s;
        if(m2 > r2) return false;

        if(l2 > r2) t = (s - sqrt(r2 - m2)) * idnorm;
        else t = (s + sqrt(r2 - m2)) * idnorm;

        if(normal) *normal = (t * dir - pos) / radius;

        return true;
    }

    inline bool sphere_raycast_full(double radius, const clam::Vec3d& pos, const clam::Vec3d& dir, double& t, clam::Vec3d* normal = nullptr){
        double s  = dot(pos, dir);
        double l2 = dot(pos, pos);
        double r2 = radius * radius;

        if(s < 0.0 && l2 > r2) return false;

        double idnorm = 1.0 / dir.length();
        s *= idnorm;
        double m2 = l2 - s * s;
        if(m2 > r2) return false;

        if(l2 > r2) t = (s - sqrt(r2 - m2)) * idnorm;
        else t = (s + sqrt(r2 - m2)) * idnorm;

        if(normal) *normal = (t * dir - pos) / radius;

        return true;
    }

    //NOTE: Assume normalized ray direction.
    //TODO: Don't assume normalized ray direction.
    inline bool ellipsoid_raycast(const clam::Vec3d& radius, const clam::Vec3d& pos, const clam::Vec3d& dir, double& t, clam::Vec3d* normal = nullptr){
        clam::Vec3d p = pos / radius;
        clam::Vec3d q = dir / radius;

        double pq = dot(p, q);
        double pp = dot(p, p);

        if(pq < 0.0 && pp > 1.0) return false;

        double qq = dot(q, q);
        double D = pq * pq - qq * (pp - 1.0);

        if(D < 0.0) return false;

        if(pp > 1.0) t = -(0.5 * pq + sqrt(D)) / qq;
        else t = -(0.5 * pq - sqrt(D)) / qq;

        if(normal){
            *normal = (t * dir - pos) / (radius * radius);
            *normal /= normal->length();
        }

        return true;
    }

    inline bool AABB_raycast(const clam::Vec3d& aabb_min, const clam::Vec3d& aabb_max, const clam::Vec3d& pos, const clam::Vec3d& dir, double &t, clam::Vec3d* normal = nullptr){
        clam::Vec3d invDir = 1.0 / dir;
        clam::Vec3d ta = (aabb_min - pos) * invDir;
        clam::Vec3d tb = (aabb_max - pos) * invDir;

        double t_min = std::max(std::max(std::min(ta[0], tb[0]), std::min(ta[1], tb[1])), std::min(ta[2], tb[2]));
        double t_max = std::min(std::min(std::max(ta[0], tb[0]), std::max(ta[1], tb[1])), std::max(ta[2], tb[2]));

        if(t_max < 0.0 || t_min > t_max) return false;

        //TODO: There must be a better way to include the calculation of the normal
        double max_dist = 0.0;
        clam::Vec3d point = pos + t * dir - 0.5 * (aabb_min + aabb_max);
        for(int i = 0; i < 3; ++i){
            if(abs(point[i]) > max_dist){
                max_dist = abs(point[i]);
                if(normal){
                    *normal = clam::Vec3d(0.0);
                    (*normal)[i] = (max_dist > 0.0)? 1.0: -1.0;
                }
            }
        }

        t = t_min;
        return true;
    }

    //NOTE: Assumes that the ray is inside the box. Returns the new cell offset.
    inline int cell_raycast(const clam::Vec3d& cellSize, const clam::Vec3d& rpos, const clam::Vec3d& dir, double& t){
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

}

#endif
