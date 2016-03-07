#ifndef __SERIALIZE_BOUNDING_VOLUME_H
#define __SERIALIZE_BOUNDING_VOLUME_H

#include "archive.h"
#include "bounding_volume_variant_fwd.h"

void serialize(Archive& ar, const bounding_volume::Variant& shape);
void deserialize(Archive& ar, bounding_volume::Variant* shape);

#endif
