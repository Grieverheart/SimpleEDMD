#include "serialization/shape_variant.h"
#include "serialization/common.h"
#include "shape/variant.h"

namespace{
    class ShapeSerializationVisitor: public boost::static_visitor<>{
    public:
        ShapeSerializationVisitor(Archive& ar):
            ar_(ar)
        {}

        template<typename T>
        void operator()(const T& shape)const{
            int shape_type = shape::index<T>::value;
            ar_.write(&shape_type, sizeof(shape_type));
            serialize(ar_, shape);
        }

    private:
        Archive& ar_;
    };
}

void serialize(Archive& ar, const shape::Variant& shape){
    boost::apply_visitor(ShapeSerializationVisitor(ar), shape);
}

void deserialize(Archive& ar, shape::Variant** shape_ptr){
    int shape_type;
    deserialize(ar, &shape_type);
    //TODO: Is there a better way to do this?
    switch(shape_type){
        case shape::index<shape::Polyhedron>::value:{
            shape::Polyhedron poly;
            deserialize(ar, &poly);
            *shape_ptr = new shape::Variant(poly);
        } break;
        case shape::index<shape::Sphere>::value:{
            shape::Sphere sph;
            deserialize(ar, &sph);
            *shape_ptr = new shape::Variant(sph);
        } break;
        case shape::index<shape::Box>::value:{
            shape::Box box;
            deserialize(ar, &box);
            *shape_ptr = new shape::Variant(box);
        } break;
        case shape::index<shape::Cone>::value:{
            shape::Cone cone;
            deserialize(ar, &cone);
            *shape_ptr = new shape::Variant(cone);
        } break;
        case shape::index<shape::Cylinder>::value:{
            shape::Cylinder cyl;
            deserialize(ar, &cyl);
            *shape_ptr = new shape::Variant(cyl);
        } break;
        //TODO: Implement panic handling.
        default: break;
    }
}
