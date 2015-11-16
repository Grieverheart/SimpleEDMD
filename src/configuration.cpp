#include "configuration.h"
#include "shape/variant.h"

Configuration::~Configuration(void){
    for(auto shape: shapes_) delete shape;
}

Configuration::Configuration(const CubicPBC& pbc, const std::vector<Particle>&  particles, const std::vector<shape::Variant*>& shapes):
    particles_(particles),
    pbc_(pbc)
{
    shapes_.resize(shapes.size());
    for(size_t i = 0; i < shapes.size(); ++i) shapes_[i] = new shape::Variant(*shapes[i]);
}

//TODO: Unsafe when config is not initialized
Configuration::Configuration(const Configuration& other):
    particles_(other.particles_),
    pbc_(other.pbc_)
{
    shapes_.resize(other.shapes_.size());
    for(size_t i = 0; i < other.shapes_.size(); ++i) shapes_[i] = new shape::Variant(*other.shapes_[i]);
}
