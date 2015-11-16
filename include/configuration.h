#ifndef __CONFIGURATION_H
#define __CONFIGURATION_H

#include "shape/variant_fwd.h"
#include "particle.h"
#include "BoundaryCondition.h"
#include <vector>

struct Configuration{
    Configuration(void){}
    Configuration(const CubicPBC& pbc, const std::vector<Particle>&  particles, const std::vector<shape::Variant*>& shapes);
    Configuration(const Configuration&);

    ~Configuration(void);

    std::vector<Particle> particles_;
    std::vector<shape::Variant*> shapes_;
    CubicPBC pbc_;
};

#endif
