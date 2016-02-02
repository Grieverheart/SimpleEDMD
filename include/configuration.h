#ifndef __CONFIGURATION_H
#define __CONFIGURATION_H

#include "shape/variant_fwd.h"
#include "serialization/archive.h"
#include "particle.h"
#include "BoundaryCondition.h"
#include <vector>

struct Configuration{
    Configuration(void){}
    Configuration(const RectangularPBC& pbc, const std::vector<Particle>&  particles, const std::vector<shape::Variant*>& shapes);
    Configuration(const Configuration&);

    ~Configuration(void);

    void serialize(Archive&)const;

    std::vector<Particle> particles_;
    std::vector<shape::Variant*> shapes_;
    RectangularPBC pbc_;
};

#endif
