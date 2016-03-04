#include "serialization/shape_variant.h"
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
