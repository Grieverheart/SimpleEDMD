#ifndef __VEC_H
#define __VEC_H

template<typename T>
struct Vec3{
    Vec3(T xx, T yy, T zz):
        x(xx), y(yy), z(zz)
    {}

    Vec3(T a):
        x(a), y(a), z(a)
    {}

    T& operator[](size_t i){
        return (&x)[i];
    }

    const T& operator[](size_t i)const{
        return (&x)[i];
    }

    Vec3<T> operator+(Vec3<T> other)const{
        return Vec3<T>(x + other.x, y + other.y, z + other.z);
    }

    Vec3<T> operator*(const T& scalar)const{
        return Vec3<T>(x * scalar, y * scalar, z * scalar);
    }

    Vec3<T>& operator+=(const Vec3<T>& other){
        return *this = *this + other;
    }

    struct{T x, y, z;};
};

typedef Vec3<double> Vec3d;
typedef Vec3<float>  Vec3f;
typedef Vec3<int>    Vec3i;

#endif
