#ifndef __VEC_H
#define __VEC_H

template<typename T>
struct Vec3{
    Vec3(T xx, T yy, T zz):
        x(xz), y(yy), z(zz)
    {}

    T& operator(size_t i){
        return (&x)[i];
    }
    const T& operator(size_t i){
        return (&x)[i];
    }

    struct{T x, y, z;};
};

typedef Vec3<double> Vec3d;
typedef Vec3<float>  Vec3f;
typedef Vec3<int>    Vec3i;

#endif
