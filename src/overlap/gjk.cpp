#include "overlap/gjk.h"
#include "shape/convex.h"
#include <cstdio>
#include "particle.h"

namespace overlap{

    namespace{
        using uchar = unsigned char;
        using uint  = unsigned int;

        using namespace clam;
        // Lookup table which tells us at which positions in the simplex array
        // to get our points a, b, c, d. I.e. if our bits are 0111 -> 7 -> {0, 1, 2}
        const uchar p_pos[16][3] =
        {
            {0, 0, 0}, {0, 0, 0}, {1, 0, 0}, {0, 1, 0},
            {2, 0, 0}, {0, 2, 0}, {1, 2, 0}, {0, 1, 2},
            {3, 0, 0}, {0, 3, 0}, {1, 3, 0}, {0, 1, 3},
            {2, 3, 0}, {0, 2, 3}, {1, 2, 3}, {0, 0, 0}
        };

        const uchar s_pos[] = {0, 0, 1, 0, 2, 0, 0, 0, 3}; //Lookup table for single enabled bit position
        //________________________^__^_____^___________^

        class Simplex{
        public:
            Simplex(void):
                max_vert2(0.0), bits_(0), last_sb_(0), size_(0)
            {}

            uchar size(void)const{
                return size_;
            }

            void add_point(const Vec3d& point){
                uchar b = ~bits_; //Flip bits
                b &= -b; //Last set (available) bit
                uchar pos = s_pos[b]; //Get the bit position from the lookup table
                last_sb_ = pos;
                bits_ |= b; //Insert the new bit
                ++size_;
                p_[pos] = point;
                double l2 = point.length2();
                if(l2 > max_vert2) max_vert2 = l2;
            }

            void remove_point(int p){
                bits_ ^= (1 << p); //Erase the bit at position p
                --size_;
            }

            const Vec3d& get_last_point(void)const{
                return p_[last_sb_];
            }

            bool contains(const Vec3d& point){
                uchar bits = bits_;
                for(int i = 0; i < 4; ++i, bits >>= 1){
                    if((bits & 1) && (p_[i] == point)) return true;
                }
                return false;
                //const uchar* pos = p_pos[(bits_ ^ (1 << last_sb_))];
                //for(int i = 0; i < size_ - 1; ++i){
                //    if(p_[pos[i]] == point) return true;
                //}
                //if(p_[last_sb_] == point) return true;
                //return false;
            }

            double max_vertex(void)const{
                return max_vert2;
            }

            void translate(const Vec3d& dr){
                //for(int k = 0; k < 4; ++k) p_[k] += dr;
                max_vert2 = 0.0;
                uchar bits = bits_;
                for(int i = 0; i < 4; ++i, bits >>= 1){
                    if(bits & 1){
                        p_[i] += dr;
                        if(p_[i].length2() > max_vert2) max_vert2 = p_[i].length2();
                    }
                }
            }

            void print(void)const{
                uchar bits = bits_;
                for(int i = 0; i < 4; ++i, bits >>= 1){
                    if(bits & 1) printf("%d: %f, %f, %f\n", i, p_[i][0], p_[i][1], p_[i][2]);
                }
            }

            void closest(Vec3d& dir);
            bool contains_origin(Vec3d& dir);

        private:
            Vec3d p_[4];//up to 4 points / 3-Simplex
            //Vec3d diffs_[4];
            double max_vert2;
            uchar bits_;
            uchar last_sb_;
            uchar size_;
        };

        inline void Simplex::closest(Vec3d& dir){
            ///////////////////////////////////////////////
            //  Check if the origin is contained in the  //
            //  Minkowski sum.                           //
            ///////////////////////////////////////////////
            switch(size_){
            case 4:
            {
                const uchar* pos = p_pos[(bits_ ^ (1 << last_sb_))];

                const Vec3d& a = p_[last_sb_];
                const Vec3d& b = p_[pos[0]];
                const Vec3d& c = p_[pos[1]];
                const Vec3d& d = p_[pos[2]];

                Vec3d ab = b - a;
                Vec3d ac = c - a;
                Vec3d ad = d - a;

                ////////////////////* Vertex Case *///////////////////

                bool after_ab = (dot(ab, a) > 0.0);
                bool after_ac = (dot(ac, a) > 0.0);
                bool after_ad = (dot(ad, a) > 0.0);

                if(after_ab && after_ac && after_ad){
                    dir = -a; //Take direction passing through origin
                    remove_point(pos[0]);
                    remove_point(pos[1]);
                    remove_point(pos[2]);
                    break;
                }

                ////////////////////* Edge Cases *///////////////////

                /* Find the triangle face normals with the correct sign (pointing outward) */
                Vec3d abxac = cross(ab, ac);
                Vec3d abxad = cross(ab, ad);

                /* ab Edge case */
                Vec3d abPerp1 = cross(abxac, ab);
                Vec3d abPerp2 = cross(abxad, ab);
                bool abPerp1Pos = (dot(abPerp1, a) > 0.0);
                bool abPerp2Pos = (dot(abPerp2, a) > 0.0);
                // The origin must be inside the space defined by the intersection
                // of two half-space normal to the adjacent faces abc, abd
                if(abPerp1Pos && abPerp2Pos && !after_ab){
                    dir = cross(cross(a, ab), ab);
                    remove_point(pos[1]);
                    remove_point(pos[2]);
                    break;
                }

                Vec3d acxad = cross(ac, ad);

                /* ac Edge case */
                Vec3d acPerp1 = cross(acxad, ac);
                Vec3d acPerp2 = cross(ac, abxac);
                bool acPerp1Pos = (dot(acPerp1, a) > 0.0);
                bool acPerp2Pos = (dot(acPerp2, a) > 0.0);
                // The origin must be inside the space defined by the intersection
                // of two half-space normal to the adjacent faces abc, acd
                if(acPerp1Pos && acPerp2Pos && !after_ac){
                    dir = cross(cross(a, ac), ac);
                    remove_point(pos[0]);
                    remove_point(pos[2]);
                    break;
                }

                /* ad Edge case */
                Vec3d adPerp1 = cross(ad, abxad);
                Vec3d adPerp2 = cross(ad, acxad);
                bool adPerp1Pos = (dot(adPerp1, a) > 0.0);
                bool adPerp2Pos = (dot(adPerp2, a) > 0.0);
                // The origin must be inside the space defined by the intersection
                // of two half-space normal to the adjacent faces acd, abd
                if(adPerp1Pos && adPerp2Pos && !after_ad){
                    dir = cross(cross(a, ad), ad);
                    remove_point(pos[0]);
                    remove_point(pos[1]);
                    break;
                }

                ////////////////////* Face Cases *///////////////////

                /* On abc side */
                Vec3d abcPerp = (dot(abxac, ad) > 0.0)? -abxac: abxac;
                // The origin should be on abc's side and between the half-spaces defined by ac and ab (normal to abc)
                if((dot(abcPerp, a) < 0.0) && !abPerp1Pos && !acPerp2Pos){
                    /* Remove point d */
                    remove_point(pos[2]);
                    dir = (dot(ad, abxac) > 0.0)? -abxac: abxac;
                    break;
                }

                /* On abd side */
                Vec3d abdPerp = (dot(abxad, ac) > 0.0)? -abxad: abxad;
                // The origin should be on abd's side and between the half-spaces defined by ab and ad (normal to abd)
                if((dot(abdPerp, a) < 0.0) && !abPerp2Pos && !adPerp1Pos){
                    /* Remove point c */
                    remove_point(pos[1]);
                    dir = (dot(ac, abxad) > 0.0)? -abxad: abxad;
                    break;
                }

                /* On acd side */
                Vec3d acdPerp = (dot(acxad, ab) > 0.0)? -acxad: acxad;
                // The origin should be on acd's side and between the half-spaces defined by ac and ad (normal to acd)
                if((dot(acdPerp, a) < 0.0) && !acPerp1Pos && !adPerp2Pos){
                    /* Remove point b */
                    remove_point(pos[0]);
                    dir = (dot(ab, acxad) > 0.0)? -acxad: acxad;
                    break;
                }

                /* 'else' should only be when the origin is inside the tetrahedron */
                break;
            }
            case 3:
            {
                const uchar* pos = p_pos[(bits_ ^ (1 << last_sb_))];

                const Vec3d& a = p_[last_sb_];
                const Vec3d& b = p_[pos[0]];
                const Vec3d& c = p_[pos[1]];

                Vec3d ab = b - a;
                Vec3d ac = c - a;

                bool after_ab = (dot(ab, a) > 0.0);
                bool after_ac = (dot(ac, a) > 0.0);

                /////////////////////* Vertex Case *///////////////////
                if(after_ab && after_ac){
                    dir = -a; //Take direction passing through origin
                    remove_point(pos[0]);
                    remove_point(pos[1]);
                    break;
                }

                ////////////////////* Edge Cases *///////////////////

                Vec3d abxac = cross(ab, ac);

                /* Origin on the outside of triangle and close to ab */
                Vec3d abPerp = cross(ab, abxac);
                if(dot(abPerp, a) < 0.0 && !after_ab){
                    dir = cross(cross(a, ab), ab);
                    /* Remove Point c */
                    remove_point(pos[1]);
                    break;
                }

                /* Origin on the outside of triangle and close to ac */
                Vec3d acPerp = cross(abxac, ac);
                if(dot(acPerp, a) < 0.0 && !after_ac){
                    dir = cross(cross(a, ac), ac);
                    /* Remove Point b */
                    remove_point(pos[0]);
                    break;
                }

                /////////////////////* Face Case *///////////////////
                dir = (dot(abxac, a) > 0.0)? -abxac: abxac;
                break;
            }
            case 2:
            {
                const uchar* pos = p_pos[(bits_ ^ (1 << last_sb_))];

                const Vec3d& a = p_[last_sb_];
                const Vec3d& b = p_[pos[0]];

                Vec3d  ab = b - a;

                if(dot(ab, a) > 0.0){
                    dir = -a; //Take direction passing through origin
                    remove_point(pos[0]);
                    break;
                }

                dir = cross(cross(a, ab), ab);
                break;
            }
            case 1:
            {
                const Vec3d& a = p_[last_sb_];
                dir = -a;
                break;
            }
            default: break;
            }
        }

        inline bool Simplex::contains_origin(Vec3d& dir){
            ///////////////////////////////////////////////
            //  Check if the origin is contained in the  //
            //  Minkowski sum.                           //
            ///////////////////////////////////////////////
            switch(size_){
            case 4:
            {
                const uchar* pos = p_pos[(bits_ ^ (1 << last_sb_))];

                const Vec3d& a = p_[last_sb_];
                const Vec3d& b = p_[pos[0]];
                const Vec3d& c = p_[pos[1]];
                const Vec3d& d = p_[pos[2]];

                Vec3d ab = b - a;
                Vec3d ac = c - a;
                Vec3d ad = d - a;

                ////////////////////* Face Cases *///////////////////

                /* On abc side */
                Vec3d abxac = cross(ab, ac);
                Vec3d abcPerp = (dot(abxac, ad) > 0.0)? -abxac: abxac;
                Vec3d abPerp1 = cross(abxac, ab);
                Vec3d acPerp2 = cross(ac, abxac);
                bool abPerp1Pos = (dot(abPerp1, a) > 0.0);
                bool acPerp2Pos = (dot(acPerp2, a) > 0.0);
                // The origin should be on abc's side and between the half-spaces defined by ac and ab (normal to abc)
                if((dot(abcPerp, a) < 0.0) && !abPerp1Pos && !acPerp2Pos){
                    /* Remove point d */
                    remove_point(pos[2]);
                    dir = abcPerp;
                    break;
                }

                /* On abd side */
                Vec3d abxad = cross(ab, ad);
                Vec3d abdPerp = (dot(abxad, ac) > 0.0)? -abxad: abxad;
                Vec3d abPerp2 = cross(abxad, ab);
                Vec3d adPerp1 = cross(ad, abxad);
                bool abPerp2Pos = (dot(abPerp2, a) > 0.0);
                bool adPerp1Pos = (dot(adPerp1, a) > 0.0);
                // The origin should be on abd's side and between the half-spaces defined by ab and ad (normal to abd)
                if((dot(abdPerp, a) < 0.0) && !abPerp2Pos && !adPerp1Pos){
                    /* Remove point c */
                    remove_point(pos[1]);
                    dir = abdPerp;
                    break;
                }

                /* On acd side */
                Vec3d acxad = cross(ac, ad);
                Vec3d acdPerp = (dot(acxad, ab) > 0.0)? -acxad: acxad;
                Vec3d acPerp1 = cross(acxad, ac);
                Vec3d adPerp2 = cross(ad, acxad);
                bool acPerp1Pos = (dot(acPerp1, a) > 0.0);
                bool adPerp2Pos = (dot(adPerp2, a) > 0.0);
                // The origin should be on acd's side and between the half-spaces defined by ac and ad (normal to acd)
                if((dot(acdPerp, a) < 0.0) && !acPerp1Pos && !adPerp2Pos){
                    /* Remove point b */
                    remove_point(pos[0]);
                    dir = acdPerp;
                    break;
                }

                ////////////////////* Edge Cases *///////////////////

                /* ab Edge case */
                // The origin must be inside the space defined by the intersection
                // of two half-space normal to the adjacent faces abc, abd
                if(abPerp1Pos && abPerp2Pos){
                    dir = cross(cross(a, ab), ab);
                    remove_point(pos[1]);
                    remove_point(pos[2]);
                    break;
                }

                /* ac Edge case */
                // The origin must be inside the space defined by the intersection
                // of two half-space normal to the adjacent faces abc, acd
                if(acPerp1Pos && acPerp2Pos){
                    dir = cross(cross(a, ac), ac);
                    remove_point(pos[0]);
                    remove_point(pos[2]);
                    break;
                }

                /* ad Edge case */
                // The origin must be inside the space defined by the intersection
                // of two half-space normal to the adjacent faces acd, abd
                if(adPerp1Pos && adPerp2Pos){
                    dir = cross(cross(a, ad), ad);
                    remove_point(pos[0]);
                    remove_point(pos[1]);
                    break;
                }

                /* 'else' should only be when the origin is inside the tetrahedron */
                return true;
            }
            case 3:
            {
                const uchar* pos = p_pos[(bits_ ^ (1 << last_sb_))];

                const Vec3d& a = p_[last_sb_];
                const Vec3d& b = p_[pos[0]];
                const Vec3d& c = p_[pos[1]];

                Vec3d ab = b - a;
                Vec3d ac = c - a;
                Vec3d abxac = cross(ab, ac);

                ////////////////////* Edge Cases *///////////////////

                /* Origin on the outside of triangle and close to ab */
                Vec3d abPerp = cross(ab, abxac);
                if(dot(abPerp, a) < 0.0){
                    dir = cross(cross(a, ab), ab);
                    /* Remove Point c */
                    remove_point(pos[1]);
                    break;
                }

                /* Origin on the outside of triangle and close to ac */
                Vec3d acPerp = cross(abxac, ac);
                if(dot(acPerp, a) < 0.0){
                    dir = cross(cross(a, ac), ac);
                    /* Remove Point b */
                    remove_point(pos[0]);
                    break;
                }

                /////////////////////* Face Case *///////////////////
                dir = (dot(abxac, a) > 0.0)? -abxac: abxac;
                break;
            }
            case 2:
            {
                const uchar* pos = p_pos[(bits_ ^ (1 << last_sb_))];

                const Vec3d& a = p_[last_sb_];
                const Vec3d& b = p_[pos[0]];

                Vec3d  ab = b - a;

                dir = cross(cross(a, ab), ab);
                break;
            }
            case 1:
            {
                const Vec3d& a = p_[last_sb_];
                dir = -a;
                break;
            }
            default: break;
            }

            return false;
        }

    };//namespace

    double gjk_distance(
        const Particle& pa, const shape::Convex& a,
        const Particle& pb, const shape::Convex& b
    ){
        auto dir = Vec3d(0.0);
        //auto dir = pb.pos - pa.pos;
        Simplex S;

        uint fail_safe = 0;

        auto inv_rot_a = pa.rot.inv();
        auto inv_rot_b = pb.rot.inv();

        while(fail_safe++ < 20){
            auto vertex_a = pa.pos + pa.rot.rotate(pa.size * a.support(inv_rot_a.rotate(dir)));
            auto vertex_b = pb.pos + pb.rot.rotate(pb.size * b.support(inv_rot_b.rotate(-dir)));
            Vec3d new_point = vertex_a - vertex_b;

            const Vec3d& last = S.get_last_point();
            //The condition is correct, it can be derived from |v|^2 - v . w < |v|^2 * eps
            //Seems like -1.0e-10 is also ok, but should make sure.
            if(S.contains(new_point) || dot(dir, new_point) - dot(dir, last) < -1.0e-10 * dot(dir, last)){
                return -dot(dir, last) / dir.length();
            }

            S.add_point(new_point);

            S.closest(dir);
            if(S.size() == 4) return 0.0;
        }

        printf("Encountered error in GJK distance: Infinite Loop.\n Direction (%f, %f, %f)\n", dir[0], dir[1], dir[2]);
        return 0.0;
    }

    bool gjk_boolean(
        const Particle& pa, const shape::Convex& a,
        const Particle& pb, const shape::Convex& b,
        double feather
    ){
        auto dir = pb.pos - pa.pos;
        Simplex S;

        uint fail_safe = 0;

        auto inv_rot_a = pa.rot.inv();
        auto inv_rot_b = pb.rot.inv();

        while(fail_safe++ < 30){
            auto vertex_a = pa.pos + pa.rot.rotate(pa.size * a.support(inv_rot_a.rotate(dir)));
            auto vertex_b = pb.pos + pb.rot.rotate(pb.size * b.support(inv_rot_b.rotate(-dir)));
            Vec3d new_point = vertex_a - vertex_b + ((feather > 0.0)? (feather / dir.length()) * dir: Vec3d(0.0));
            double dn = dot(dir, new_point);
            if(dn < 0.0 || S.contains(new_point)) return false;
            S.add_point(new_point);
            if(S.contains_origin(dir)) return true;
        }

        printf("Encountered error in GJK boolean: Infinite Loop.\n Direction (%f, %f, %f)\n", dir[0], dir[1], dir[2]);

        return true;
    }

    bool gjk_raycast(
        const Particle& pa, const shape::Convex& a,
        const Particle& pb, const shape::Convex& b,
        const Vec3d& ray_dir, double& distance,
        clam::Vec3d& normal
    )
    {
        Vec3d dir = pb.pos - pa.pos;
        Simplex S;

        uint fail_safe = 0;

        Vec3d x(0.0);
        double lambda = 0.0;

        auto inv_rot_a = pa.rot.inv();
        auto inv_rot_b = pb.rot.inv();

        while(fail_safe++ < 100){
            auto vertex_a = pa.pos + pa.rot.rotate(pa.size * a.support(inv_rot_a.rotate(dir)));
            auto vertex_b = pb.pos + pb.rot.rotate(pb.size * b.support(inv_rot_b.rotate(-dir)));
            Vec3d new_point = vertex_a - vertex_b;

            Vec3d new_point_trans = new_point - x;

            if(dot(dir, new_point_trans) < 0.0){
                if(dot(dir, ray_dir) <= 0.0) return false;

                double delta = dot(dir, new_point_trans) / dot(dir, ray_dir);
                lambda -= delta;
                if(lambda > distance) return false;
                x = -lambda * ray_dir;
                S.translate(-delta * ray_dir);
            }

            //if(sqr(dot(dir, new_point_trans)) < 2.41e-30 * S.max_vertex() * dot(dir, dir)){
            //    distance = lambda;
            //    return true;
            //}

            S.add_point(new_point_trans);
            S.closest(dir);

            if(S.size() == 4){
                distance = lambda;
                normal = dir / dir.length();
                return true;
            }
        }
        distance = lambda;
        normal = dir / dir.length();
        ///printf("Encountered error in GJK raycast: Infinite Loop.\n Direction (%f, %f, %f)\n", dir[0], dir[1], dir[2]);
        //auto test = pb.pos - pa.pos;
        //for(int i = 0; i < n_ls; ++i) printf("%f, %f, %f\n", lambdas[i][0], lambdas[i][1], lambdas[i][2]);
        //printf("Encountered error in GJK raycast: Infinite Loop.\n Direction (%f, %f, %f)\n", test[0], test[1], test[2]);
        //printf("%f: %f, %f, %f\n", distance, ray_dir[0], ray_dir[1], ray_dir[2]);
        return true;
    }
}
