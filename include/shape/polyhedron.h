#ifndef EDMD_SHAPE_POLYHEDRON_H
#define EDMD_SHAPE_POLYHEDRON_H

#include "convex.h"
#include <cstddef>
#include <vector>

namespace shape{

    //TODO: In the future, we might allow a size parameter for the shape itself, then
    //we can maybe make the Particle.size member a pointer to double so that we can point
    //to the shape's size instead
    //TODO: Get rid of unsigned int. Better use something like uint_fast32_t
    class Polyhedron: public Convex{
    public:
        explicit Polyhedron(const std::vector<clam::Vec3d>& vertices, const std::vector<std::vector<unsigned int>>& faces, const char* source = nullptr);
        Polyhedron(const Polyhedron&);
        ~Polyhedron(void);
        //TODO: Add void set_source(const char*) and const char* get_source(void) functions
        //so that we can handle the storage of an optional filename from which the shape
        //was loaded from.
        double in_radius(void)const{
            return in_radius_;
        }
        double out_radius(void)const{
            return out_radius_;
        }
        double volume(void)const{
            return volume_;
        }
        const std::vector<clam::Vec3d>& vertices(void)const{
            return vertices_;
        }
        const std::vector<std::vector<unsigned int>>& faces(void)const{
            return faces_;
        }
        const clam::Vec3d& get_vertex(size_t idx)const{
            return vertices_[idx];
        }
        const std::vector<unsigned int>& get_vertex_nbs(size_t idx)const{
            return vert_neighbors[idx];
        }

        const char* get_source(void)const{
            return source_;
        }

        clam::Vec3d support(const clam::Vec3d&)const;
        double max_vert_dist2(const clam::Vec3d& pos, const clam::Quatd& rot)const;
    private:
        //@note: perhaps hold the square of these values
        char* source_;
        double in_radius_;
        double out_radius_;
        double volume_;
        std::vector<clam::Vec3d> vertices_;
        std::vector<std::vector<unsigned int>> vert_neighbors;
        std::vector<std::vector<unsigned int>> faces_;
    };

    //@note: As a future reminder, we would like to also implement a shape union.
    //class UnionShape{
    //    double out_radius(void)const;
    //    std::vector<Shape> shapes_;
    //};

    inline clam::Vec3d Polyhedron::support(const clam::Vec3d& dir)const{
        using idx_t = std::vector<unsigned int>::size_type;
        idx_t curr = 0;
        double p = 0.0;
        double max = clam::dot(vertices_[0], dir);
        for(idx_t i = 1; i < vertices_.size(); ++i){
            p = clam::dot(vertices_[i], dir);
            if(p > max){
                curr = i;
                max = p;
            }
        }
        return vertices_[curr];
        //using idx_t = std::vector<unsigned int>::size_type;
        //unsigned int next = 0, last = 0, curr = 0;
        //idx_t curr = 0;
        //double p = 0.0;
        //double max = clam::dot(vertices_[0], dir);
        //for(;;){
        //    for(idx_t vid = 0; vid < vert_neighbors[curr].size(); ++vid){
        //        next = vert_neighbors[curr][vid];
        //        if(next != last){
        //            p = clam::dot(vertices_[next], dir);
        //            if(p > max){
        //                max = p;
        //                last = curr;
        //                curr = next;
        //                break;
        //            }
        //        }
        //        if(vid == vert_neighbors[curr].size() - 1) return vertices_[curr];
        //    }
        //}
    }

    inline double Polyhedron::max_vert_dist2(const clam::Vec3d& pos, const clam::Quatd& rot)const{
        using idx_t = std::vector<unsigned int>::size_type;
        double p = 0.0;
        double max = (rot.rotate(vertices_[0]) + pos).length2();
        for(idx_t i = 1; i < vertices_.size(); ++i){
            p = (rot.rotate(vertices_[i]) + pos).length2();
            if(p > max) max = p;
        }
        return max;
        //using idx_t = std::vector<unsigned int>::size_type;
        //unsigned int next = 0, last = 0, curr = 0;
        //double p = 0.0;
        //double max = (rot.rotate(vertices_[0]) + pos).length2();
        //for(;;){
        //    for(idx_t vid = 0; vid < vert_neighbors[curr].size(); ++vid){
        //        next = vert_neighbors[curr][vid];
        //        if(next != last){
        //            p = (rot.rotate(vertices_[next]) + pos).length2();
        //            if(p > max){
        //                max = p;
        //                last = curr;
        //                curr = next;
        //                break;
        //            }
        //        }
        //        if(vid == vert_neighbors[curr].size() - 1) return max;
        //    }
        //}
    }

}
#endif
