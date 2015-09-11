#ifndef EDMD_OBJ_LOADER_H
#define EDMD_OBJ_LOADER_H

#include <vector>
#include "clam.h"

bool load_obj(const char* filepath, std::vector<clam::Vec3d>& vertices, std::vector<std::vector<unsigned int>>& faces);

#endif
