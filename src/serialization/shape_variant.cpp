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
            int shape_idx = shape::index<T>::value;
            ar_.write(&shape_idx, sizeof(shape_idx));
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
    int shape_idx;
    deserialize(ar, &shape_idx);
    //TODO: Is there a better way to do this?
    //NOTE: We can extract a type from the index; shape::variant_index_type<index>::type.
    //Then,
    //using shape_type = shape::variant_index_type<shape_idx>::type;
    //shape_type shape;
    //deserialize(ar, &shape);
    //*shape_ptr = new shape::Variant(shape);
    switch(shape_idx){
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
        case shape::index<shape::LeafCylinder>::value:{
            shape::LeafCylinder leaf;
            deserialize(ar, &leaf);
            *shape_ptr = new shape::Variant(leaf);
        } break;
        //TODO: Implement panic handling.
        default: break;
    }
}
