#ifndef CLAM_H
#define CLAM_H

#ifdef __cplusplus

#include "serialization/archive.h"
#include <cmath>

namespace clam{
    namespace{

        template<typename T>
        static inline T sqr(T val){
            return val * val;
        }

        // Fast and accurate cube-root by Kahan
        inline double fast_cbrt(double k){
            const unsigned int B1 = 715094163u;
            double t = 0.0;
            unsigned int* pt = (unsigned int*) &t;
            unsigned int* px = (unsigned int*) &k;
            pt[1] = px[1] / 3 + B1;
            double x = t;
            for(int i = 0; i < 2; i++){
                double tri = x*x*x;
                x = x * ((tri + k + k) / (tri + tri + k));
            }
            return x;
        }
        
        inline float fast_cbrt(float k){
            float f = k;
            unsigned int* p = (unsigned int*) &f;
            *p = *p / 3 + 709921077u;
            return f;
        }
    };
	/* Vector */
    template<typename T>
	class Vec3{
	private:
		T data_[3];
	public:
		constexpr Vec3(void): data_{0}
        {}

		constexpr Vec3(T a):
            data_{a, a, a}
        {}

		constexpr Vec3(T x, T y, T z):
            data_{x, y, z}
        {}

		constexpr Vec3(const T* a):
            data_{a[0], a[1], a[2]}
        {}

        constexpr const T* data(void)const{
            return data_;
        }

		constexpr T operator[](unsigned int i)const{
			return data_[i];
		}

		T& operator[](unsigned int i){
			return data_[i];
		}

        constexpr bool operator==(const Vec3<T>& other)const{
            return ((other[0] == data_[0]) &&
                    (other[1] == data_[1]) &&
                    (other[2] == data_[2]));
        }

        constexpr Vec3<T> operator*(T a)const{
            return Vec3<T>(data_[0] * a, data_[1] * a, data_[2] * a);
        }

        Vec3<T>& operator*=(T a){
            data_[0] *= a;
            data_[1] *= a;
            data_[2] *= a;
            return *this;
        }

        constexpr Vec3<T> operator*(const Vec3<T>& a)const{
            return Vec3<T>(data_[0] * a[0], data_[1] * a[1], data_[2] * a[2]);
        }

        Vec3<T> operator*=(const Vec3<T>& a){
            data_[0] *= a.data_[0];
            data_[1] *= a.data_[1];
            data_[2] *= a.data_[2];
            return *this;
        }

        constexpr Vec3<T> operator/(T a)const{
            return Vec3<T>(data_[0] / a, data_[1] / a, data_[2] / a);
        }

        Vec3<T> operator/=(T a){
            data_[0] /= a;
            data_[1] /= a;
            data_[2] /= a;
            return *this;
        }

        constexpr Vec3<T> operator/(const Vec3<T>& a)const{
            return Vec3<T>(data_[0] / a[0], data_[1] / a[1], data_[2] / a[2]);
        }

        Vec3<T>& operator/=(const Vec3<T>& a){
            data_[0] /= a[0];
            data_[1] /= a[1];
            data_[2] /= a[2];
            return *this;
        }

        constexpr Vec3<T> operator+(const Vec3<T>& other)const{
            return Vec3<T>(
                data_[0] + other[0],
                data_[1] + other[1],
                data_[2] + other[2]
            );
        }

        Vec3<T>& operator+=(const Vec3<T>& other){
            data_[0] += other[0];
            data_[1] += other[1];
            data_[2] += other[2];
            return *this;
        }

        constexpr Vec3<T> operator-(const Vec3<T>& other)const{
            return Vec3<T>(
                data_[0] - other[0],
                data_[1] - other[1],
                data_[2] - other[2]
            );
        }

        constexpr Vec3<T> operator-(void)const{
            return Vec3<T>(-data_[0], -data_[1], -data_[2]);
        }

        Vec3<T>& operator-=(const Vec3<T>& other){
            data_[0] -= other[0];
            data_[1] -= other[1];
            data_[2] -= other[2];
            return *this;
        }

        constexpr T length2(void)const{
            return data_[0] * data_[0] + data_[1] * data_[1] + data_[2] * data_[2];
        }

        constexpr T length(void)const{
            return sqrt(length2());
        }

	};

    template<typename T>
    constexpr Vec3<T> cross(const Vec3<T>& a, const Vec3<T>& b){
		return Vec3<T>(
            a[1] * b[2] - a[2] * b[1],
            a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0]
        );
    }

    template<typename T>
    constexpr Vec3<T> triple(const Vec3<T>& a, const Vec3<T>& b, const Vec3<T>& c){
		return Vec3<T>(
            b[0] * dot(a, c) - a[0] * dot(b, c),
            b[1] * dot(a, c) - a[1] * dot(b, c),
            b[2] * dot(a, c) - a[2] * dot(b, c)
        );
    }

    template<typename T>
    constexpr T dot(const Vec3<T>& a, const Vec3<T>& b){
		return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
    }

    template<typename T>
    constexpr Vec3<T> operator*(T a, const Vec3<T>& b){
        return Vec3<T>(a * b[0], a * b[1], a * b[2]);
    }

    template<typename T>
    constexpr Vec3<T> operator/(T a, const Vec3<T>& b){
        return Vec3<T>(a / b[0], a / b[1], a / b[2]);
    }

	/* Quaternion */
    //TODO: Maybe it's better to reverse the order of v_, w_
    template<typename T>
	class Quat{
    private:
        Vec3<T> v_;
        T w_;
	public:
		constexpr Quat(void):
            v_(0.0), w_(0.0)
        {}

		constexpr Quat(T x, T y, T z, T w):
            v_(x, y, z), w_(w)
        {}

		constexpr Quat(const Vec3<T>& v, T w):
            v_(v), w_(w)
        {}

        constexpr Quat<T> conj(void)const{
            return Quat<T>(-v_, w_);
        }

        constexpr Vec3<T> imag(void)const{
            return v_;
        }

        constexpr T real(void)const{
            return w_;
        }

		constexpr T operator[](unsigned int i)const{
			return (i < 3)? v_[i]: w_;
		}

		T& operator[](unsigned int i){
			return (i < 3)? v_[i]: w_;
		}

        constexpr Quat<T> operator*(const Quat<T>& other)const{
            return Quat<T>(
                w_ * other.v_ + other.w_ * v_ + cross(v_, other.v_),
                w_ * other.w_ - dot(v_, other.v_)
            );
        }

        constexpr Quat<T> operator*(T a)const{
            return Quat<T>(a * v_[0], a * v_[1], a * v_[2], a * w_);
        }

        Quat<T> operator*=(T a){
            v_[0] *= a;
            v_[1] *= a;
            v_[2] *= a;
            w_ *= a;
            return *this;
        }

        constexpr Quat<T> operator/(T a)const{
            return Quat<T>(v_[0] / a, v_[1] / a, v_[2] / a, w_ / a);
        }

        Quat<T> operator/=(T a){
            v_[0] /= a;
            v_[1] /= a;
            v_[2] /= a;
            w_ /= a;
            return *this;
        }

        constexpr Quat<T> operator+(const Quat<T>& other)const{
            return Quat<T>(
                v_[0] + other.v_[0],
                v_[1] + other.v_[1],
                v_[2] + other.v_[2],
                w_ + other.w_
            );
        }

        Quat<T> operator+=(const Quat<T>& other){
            v_[0] += other.v_[0];
            v_[1] += other.v_[1];
            v_[2] += other.v_[2];
            w_ += other.w_;
            return *this;
        }

        constexpr Quat<T> operator-(const Quat<T>& other)const{
            return Quat<T>(
                v_[0] - other.v_[0],
                v_[1] - other.v_[1],
                v_[2] - other.v_[2],
                w_ - other.w_
            );
        }

        Quat<T> operator-=(const Quat<T>& other){
            v_[0] -= other.v_[0];
            v_[1] -= other.v_[1];
            v_[2] -= other.v_[2];
            w_ -= other.w_;
            return *this;
        }

        //TODO: Check if we need to divide by lengh2 if already normalized
        constexpr Quat<T> inv(void)const{
            return conj() / length2();
        }

        constexpr T length2(void)const{
            return v_[0] * v_[0] + v_[1] * v_[1] + v_[2] * v_[2] + w_ * w_;
        }

        constexpr T length(void)const{
            return sqrt(length2());
        }

        constexpr Vec3<T> rotate(const Vec3<T>& vec)const{
            return vec + 2.0 * cross(v_, cross(v_, vec) + w_ * vec);
        }

        void to_matrix(T* data)const{
            double q00 = v_[0] * v_[0];
            double q11 = v_[1] * v_[1];
            double q22 = v_[2] * v_[2];
            double q01 = v_[0] * v_[1];
            double q02 = v_[0] * v_[2];
            double q12 = v_[0] * v_[2];
            double q0w = v_[0] * w_;
            double q1w = v_[1] * w_;
            double q2w = v_[2] * w_;

            data[0] = 1.0 - 2.0 * (q11 + q22);
            data[1] = 2.0 * (q01 + q2w);
            data[2] = 2.0 * (q02 - q1w);
            data[3] = 2.0 * (q01 - q2w);
            data[4] = 1.0 - 2.0 * (q00 + q22);
            data[5] = 2.0 * (q12 + q0w);
            data[6] = 2.0 * (q02 + q1w);
            data[7] = 2.0 * (q12 - q0w);
            data[8] = 1.0 - 2.0 * (q00 + q11);
        }

        void toAxisAngle(T& angle, Vec3<T>& axis)const{
            if(w_ * w_ < 1.0){
                T s = 1.0 / sqrt(1.0 - w_ * w_);
                angle = 2.0 * acos(w_);
                axis = s * v_;
                axis /= axis.length();
            }
            else{
                angle = 0.0;
                axis = Vec3<T>(0.0, 1.0, 0.0);
            }
        }
    };

    template<typename T>
    Quat<T> fromAxisAngle(T angle, const clam::Vec3<T>& axis){
        return Quat<T>(sin(0.5 * angle) * axis, cos(0.5 * angle));
    }

    template<typename T, typename F = T(void)>
    Vec3<T> uniform_vector(F gen_rand){
        T x1, x2;
        do{
            x1 = 2.0 * gen_rand() - 1.0;
            x2 = 2.0 * gen_rand() - 1.0;
        }while(x1 * x1 + x2 * x2 >= 1.0);

        return clam::Vec3<T>(
            2.0 * x1 * sqrt(1.0 - x1 * x1 - x2 * x2),
            2.0 * x2 * sqrt(1.0 - x1 * x1 - x2 * x2),
            1.0 - 2.0 * (x1 * x1 + x2 * x2)
        );
    }

    template<typename T, typename F = T(void)>
    Quat<T> uniform_quaternion(T max_angle, F gen_rand){
        //Generate Random Angle
        T angle = max_angle * max_angle * max_angle * gen_rand();
        angle = fast_cbrt(angle);
        while(gen_rand() > sqr(sin(angle) / angle)){
            angle = max_angle * max_angle * max_angle * gen_rand();
            angle = fast_cbrt(angle);
        }
        T sina = sin(angle * 0.5);

        //TODO: Use fromAxisAngle.
        return Quat<T>(
            uniform_vector<T, F>(gen_rand) * sina,
            cos(0.5 * angle)
        );
    }

    using Vec3d = Vec3<double>;
    using Vec3f = Vec3<float>;
    using Quatd = Quat<double>;
    using Quatf = Quat<float>;
};

#endif
#endif
