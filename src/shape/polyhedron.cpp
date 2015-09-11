#include "shape/polyhedron.h"
#include <cstring>

namespace shape{

    typedef unsigned int uint;

    Polyhedron::Polyhedron(const std::vector<clam::Vec3d>& vertices, const std::vector<std::vector<uint>>& faces, const char* source):
        source_(nullptr), volume_(0.0), vertices_(vertices), faces_(faces)
    {
        if(source){
            source_ = new char[strlen(source) + 1];
            strcpy(source_, source);
        }

        using idx_t = std::vector<std::vector<uint>>::size_type;
        /* Calculate Properties */

        /* Find Edges */
        //Iterate over each face and for each next face in the list, check if they
        //share two vertices, this defines an edge.
        std::vector< std::vector<uint> > edges;
        for(idx_t fi = 0; fi < faces.size(); ++fi){

            auto normal = clam::cross(
                vertices_[faces[fi][1]] - vertices_[faces[fi][0]],
                vertices_[faces[fi][faces[fi].size() - 1]] - vertices_[faces[fi][0]]
            );
            normal /= normal.length();

            auto barycenter = clam::Vec3d(0.0);

            double area = 0.0;
            for(idx_t fvidx = 0; fvidx < faces[fi].size(); ++fvidx){
                const auto& v1 = vertices_[faces[fi][fvidx]];
                const auto& v2 = vertices_[faces[fi][(fvidx + 1) % faces[fi].size()]];
                area += 0.5 * clam::dot(normal, clam::cross(v1, v2));
                barycenter += vertices_[faces[fi][fvidx]];
            }
            volume_ += fabs(clam::dot(barycenter, normal) * area) / (3.0 * faces[fi].size());
            for(idx_t fj = fi + 1; fj < faces.size(); ++fj){
                uint fcount = 0;
                std::vector<uint> edge;
                for(auto face1_vidx: faces[fi]){
                    for(auto face2_vidx: faces[fj]){
                        if(face1_vidx == face2_vidx){
                            edge.push_back(face1_vidx);
                            ++fcount;
                        }
                    }
                    if(fcount == 2){
                        edges.push_back(edge);
                        fcount = 0;
                        edge.clear();
                    }
                }
            }
        }
        //printf("Found %lu edges\n", edges.size());

        /* Find Vertex Neighbours */
        //For all vertices, check if two edges share this vertex. If they do and it
        //isn't vertex 0, append the other vertices of these edge to the neighbor list
        for(uint vi = 0; vi < vertices_.size(); ++vi){
            std::vector<uint> neighbors;
            for(uint ei = 0; ei < edges.size(); ++ei){
                for(uint i = 0; i < 2; ++i){
                    if(edges[ei][i] == vi && edges[ei][(i + 1) % 2] != 0) neighbors.push_back(edges[ei][(i + 1) % 2]);
                }
            }
            if(!neighbors.empty()) vert_neighbors.push_back(neighbors);
        }

        /* Find the inscribed sphere radius */
        //For each face, calculate its distance from the particle's center and find the min
        double minDistance = 100.0;
        for(auto faceItr = faces.begin(); faceItr < faces.end(); faceItr++){
            clam::Vec3d p(vertices_[(*faceItr)[0]]);

            clam::Vec3d a(
                vertices_[(*faceItr)[1]][0] - vertices_[(*faceItr)[0]][0],
                vertices_[(*faceItr)[1]][1] - vertices_[(*faceItr)[0]][1],
                vertices_[(*faceItr)[1]][2] - vertices_[(*faceItr)[0]][2]
            );

            clam::Vec3d b(
                vertices_[(*faceItr)[2]][0] - vertices_[(*faceItr)[0]][0],
                vertices_[(*faceItr)[2]][1] - vertices_[(*faceItr)[0]][1],
                vertices_[(*faceItr)[2]][2] - vertices_[(*faceItr)[0]][2]
            );

            clam::Vec3d normal = clam::cross(a, b);
            double length = normal.length();
            for(int i = 0; i < 3; ++i) normal[i] /= length;
            double faceDistance = fabs(clam::dot(normal, p));

            if(faceDistance < minDistance) minDistance = faceDistance;
        }
        in_radius_ = minDistance;
        //printf("Inscribed Radius: %f\n", in_radius_);

        /* Find the circumscribed sphere radius */
        //It's just the farthest vertex from the particle's center
        double maxDistance = 0.0;
        for(auto vItr = vertices_.begin(); vItr < vertices_.end(); vItr++){
            double vertexLength = vItr->length();
            if( vertexLength > maxDistance) maxDistance = vertexLength;
        }
        out_radius_ = maxDistance;
        //printf("Circumscribed Radius: %f\n", out_radius_);
    }

    Polyhedron::Polyhedron(const Polyhedron& other):
        source_(nullptr), 
        in_radius_(other.in_radius_), out_radius_(other.out_radius_),
        volume_(other.volume_), 
        vertices_(other.vertices_), 
        vert_neighbors(other.vert_neighbors),
        faces_(other.faces_)
    {
        if(other.source_){
            source_ = new char[strlen(other.source_) + 1];
            strcpy(source_, other.source_);
        }
    }


    Polyhedron::~Polyhedron(void){
        delete[] source_;
    }

}
