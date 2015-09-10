#ifndef EDMD_SHAPE_VARIANT_H
#define EDMD_SHAPE_VARIANT_H

#include "shape/shapes.h"
#include "shape/variant_fwd.h"
#include <boost/variant.hpp>

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
