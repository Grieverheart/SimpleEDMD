#include "serialization/shape_variant.h"
#include "serialization/common.h"
#include "shape/variant.h"

namespace{
    class ShapeSerializationVisitor: public boost::static_visitor<>{
    public:
        ShapeSerializationVisitor(Archive& ar):
            ar_(ar)
        {}

        void operator()(const shape::Polyhedron& poly)const{
            int shape_type = shape::POLYHEDRON;
            ar_.write(&shape_type, sizeof(shape_type));
            serialize(ar_, poly);
        }

        void operator()(const shape::Sphere& sph)const{
            int shape_type = shape::SPHERE;
            ar_.write(&shape_type, sizeof(shape_type));
            serialize(ar_, sph);
        }

    private:
        Archive& ar_;
    };
}

void serialize(Archive& ar, const shape::Variant& shape){
    boost::apply_visitor(ShapeSerializationVisitor(ar), shape);
}

void deserialize(Archive& ar, shape::Variant* shape_ptr){
    int shape_type;
    deserialize(ar, &shape_type);
    switch(shape_type){
        case shape::POLYHEDRON:{
            shape::Polyhedron poly;
            deserialize(ar, &poly);
            shape_ptr = new shape::Variant(poly);
        } break;
        case shape::SPHERE:{
            shape::Sphere sph;
            deserialize(ar, &sph);
            shape_ptr = new shape::Variant(sph);
        } break;
        //TODO: Implement panic handling.
        default: break;
    }
}
