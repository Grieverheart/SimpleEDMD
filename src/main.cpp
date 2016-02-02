#include <cstdio>
#include <cstring>
#include <unistd.h>
#include "shape/variant.h"
#include "Simulation.h"
#include "obj_loader.h"
#include "io/config_xml.h"
#include "overlap/overlap.h"
#include "serialization/archive.h"

inline void stream_position(Particle& particle, double time){
    particle.pos += particle.vel * (time - particle.time);
}

inline void stream_rotation(Particle& particle, double time){
    clam::Vec3d ha = ((time - particle.time) * 0.5) * particle.ang_vel;
    double l = ha.length(); // magnitude
    if(l > 0.0){
        double sl, cl;
        sincos(l, &sl, &cl);
        particle.rot = clam::Quatd(ha * (sl / l), cl) * particle.rot;
    }
}

inline void update_particle(Particle& particle, double time, const RectangularPBC& pbc){
    if(particle.time < time){
        stream_position(particle, time);
        stream_rotation(particle, time);
        particle.pos = pbc.apply(particle.pos);
        particle.time = time;
    }
}

double packing_fraction(const Configuration& config){
    auto box_size = config.pbc_.getSize();
    double volume = box_size[0] * box_size[1] * box_size[2];

    //NOTE: Assume particle unit volume!
    return config.particles_.size() / volume;
}

//TODO: Add function to reset simulation statistics (av_momentum_transfer).
int main(int argc, char *argv[]){

    Simulation* sim;

    {
        Configuration config;
        xml_load_config(argv[1], config);

        //Archive ar;
        //config.serialize(ar);
        //FILE* fp = fopen("Data/config.core", "wb");
        //fwrite(ar.data(), 1, ar.size(), fp);
        //fclose(fp);

        //Scale box and positions for reaching target_pf
        {
            double target_pf = atof(argv[2]);
            double pf = packing_fraction(config);
            double scale = pow(pf / target_pf, 1.0 / 3.0);
            config.pbc_.setSize(scale * config.pbc_.getSize());
            for(auto& particle: config.particles_) particle.pos = scale * particle.pos;
        }

        sim = new Simulation(config);
    }

    PeriodicCallback output(0.01);
    output.setNextFunction([](double time){
        return time + 1.0;
    });
    output.setCallback([sim](double time){
        static int nFiles = 0;
        Configuration config = sim->get_configuration();

        for(auto& particle: config.particles_){
            update_particle(particle, time, config.pbc_);
        }

        char buff[64];
        sprintf(buff, "Data/pid%u.%06u.xml", getpid(), ++nFiles);
        xml_save_config(buff, config);
        clam::Vec3d box_size = config.pbc_.getSize();
        double volume = box_size[0] * box_size[1] * box_size[2];
        double kT = 2.0 * sim->get_average_kinetic_energy() / (3.0 * config.particles_.size());
        double pressure = (config.particles_.size() - sim->get_average_stress() / (3.0 * kT)) / volume;
        printf("%e: %f\t%f\n", time, pressure, kT);
        sim->reset_statistics();

#ifndef NDEBUG
        //Check for overlaps
        bool overlaps = false;
        for(size_t i = 0; i < config.particles_.size(); ++i){
            for(size_t j = i + 1; j < config.particles_.size(); ++j){
                Particle pa = config.particles_[i];
                Particle pb = config.particles_[j];
                pb.pos = config.pbc_.minImage(pb.pos - pa.pos);
                pa.pos = 0.0;
                if(overlap::shape_overlap(pa, *config.shapes_[pa.shape_id], pb, *config.shapes_[pb.shape_id])){
                    printf("%lu, %lu\n", i, j);
                    overlaps = true;
                }
            }
        }
        assert(overlaps == false);
#endif
    });

    sim->run(1000.0, output);

    return 0;
}
