#ifndef EDMD_SHAPE_VARIANT_H
#define EDMD_SHAPE_VARIANT_H

#include "shape/shapes.h"
#include "shape/variant_fwd.h"
#include <boost/variant.hpp>

//NOTE: Boost Variant does not implement move semantics.

namespace shape{
    namespace detail{
        template<int idx, typename T, typename... ArgsT>
        struct type_index_helper;

        template<int idx, typename T>
        struct type_index_helper<idx, T>{
            static constexpr int value = -1;
        };

        template<int idx, typename T, typename A, typename... ArgsT>
        struct type_index_helper<idx, T, A, ArgsT...>{
            static constexpr int value = (std::is_same<T, A>::value)? idx: type_index_helper<idx + 1, T, ArgsT...>::value;
        };

        template<typename T, typename... ArgsT>
        struct type_index{
            static constexpr int value = type_index_helper<0, T, ArgsT...>::value;
        };
    }

    template<typename T>
    struct index{
        static constexpr int value = detail::type_index<T, Polyhedron, Sphere, Box, Cone, Cylinder>::value;
    };
}

class ShapeOutRadiusVisitor: public boost::static_visitor<double>{
public:
    template<typename T>
    double operator()(const T& shape)const{
        return shape.out_radius();
    }

    double operator()(const shape::Sphere& sph)const{
        return sph.radius();
    }
};

inline double shape_outradius(const shape::Variant& shape){
    return boost::apply_visitor(ShapeOutRadiusVisitor(), shape);
}

//TODO: Perhaps do the following to avoid calling boost explicitly
//
//namespace shape{
//    template<class T>
//    using static_visitor = boost::static_visitor<T>;
//
//    template <typename... Args>
//    auto apply_visitor(Args&&... args) -> decltype(boost::apply_visitor(std::forward<Args>(args)...)){
//        return boost::apply_visitor(std::forward<Args>(args)...);
//    }
//}
//
//TODO: We can apply a visitor in plain C by defining a generic visitor which takes
//n_shapes function pointers i.e.
//namespace shape{
//    class ShapeGenericVisitor: public boost::static_visitor<void*>{
//    public:
//        using poly_func = void* (*)(const Polyhedron*);
//        using sph_func  = void* (*)(const Sphere*);
//
//        ShapeGenericVisitor(poly_func pf, sph_func sf):
//            pf_(pf), sf_(sf)
//        {}
//
//        void* operator()(const Polyhedron& poly)const{
//            return pf_(&poly);
//        }
//
//        void* operator()(const Sphere& sph)const{
//            return sf_(&sph);
//        }
//
//    private:
//        poly_func pf_;
//        sph_func sf_;
//    };
//}
//
//void* shape_visit(shape::Variant var, poly_func pf, sph_func sf){
//    auto visitor = ShapeGenericVisitor(pf, sf);
//    return boost::apply_visitor(visitor, var);
//}

#endif
