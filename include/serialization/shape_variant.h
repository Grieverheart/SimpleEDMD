#ifndef __SERIALIZE_SHAPE_H
#define __SERIALIZE_SHAPE_H

#include "archive.h"
#include "shape/variant_fwd.h"

void serialize(Archive& ar, const shape::Variant& shape);

#endif
