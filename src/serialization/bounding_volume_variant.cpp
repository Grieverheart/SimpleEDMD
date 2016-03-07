#include "serialization/bounding_volume_variant.h"
#include "serialization/common.h"
#include "bounding_volume_variant.h"

namespace{
    class BoundingVolumeSerializationVisitor: public boost::static_visitor<>{
    public:
        BoundingVolumeSerializationVisitor(Archive& ar):
            ar_(ar)
        {}

        void operator()(const shape::Box& box)const{
            int shape_type = bounding_volume::BOX;
            ar_.write(&shape_type, sizeof(shape_type));
            serialize(ar_, box);
        }

        void operator()(const shape::Sphere& sph)const{
            int shape_type = bounding_volume::SPHERE;
            ar_.write(&shape_type, sizeof(shape_type));
            serialize(ar_, sph);
        }

    private:
        Archive& ar_;
    };
}

void serialize(Archive& ar, const bounding_volume::Variant& shape){
    boost::apply_visitor(BoundingVolumeSerializationVisitor(ar), shape);
}

void deserialize(Archive& ar, bounding_volume::Variant* shape_ptr){
    int shape_type;
    deserialize(ar, &shape_type);
    switch(shape_type){
        case bounding_volume::BOX:{
            shape::Box box;
            deserialize(ar, &box);
            shape_ptr = new bounding_volume::Variant(box);
        } break;
        case bounding_volume::SPHERE:{
            shape::Sphere sph;
            deserialize(ar, &sph);
            shape_ptr = new bounding_volume::Variant(sph);
        } break;
        //TODO: Implement panic handling.
        default: break;
    }
}
