#ifndef EDMD_SHAPE_COMPLEX_H
#define EDMD_SHAPE_COMPLEX_H

#include "variant_fwd.h"
#include "clam.h"
#include "transformation.h"
#include <vector>

namespace shape{
    class Complex{
    public:
        struct SubShape{
            Transformation xform_;
            Variant* shape_;
        };

        Complex(void):
            out_radius_(0.0)
        {}

        template<class... ArgsT>
        void emplace_shape(ArgsT&&... args){
            shapes_.emplace_back(std::forward<ArgsT>(args)...);
            update_outradius();
        }

        const std::vector<SubShape>& shapes(void){
            return shapes_;
        }

        double out_radius(void)const{
            return out_radius_;
        }

    private:
        void update_outradius(void);

        double out_radius_;
        std::vector<SubShape> shapes_;
    };
}

#endif
