#include <cstdio>
#include <cstring>
#include <unistd.h>
#include "shape/variant.h"
#include "Simulation.h"
#include "obj_loader.h"
#include "io/config_xml.h"

inline void stream_position(Particle& particle, double time){
    particle.xform_.pos_ += particle.vel * (time - particle.time);
}

inline void stream_rotation(Particle& particle, double time){
    clam::Vec3d ha = ((time - particle.time) * 0.5) * particle.ang_vel;
    double l = ha.length(); // magnitude
    if(l > 0.0){
        double sl, cl;
        sincos(l, &sl, &cl);
        particle.xform_.rot_ = clam::Quatd(ha * (sl / l), cl) * particle.xform_.rot_;
    }
}

inline void update_particle(Particle& particle, double time, const RectangularPBC& pbc){
    if(particle.time < time){
        stream_position(particle, time);
        stream_rotation(particle, time);
        particle.xform_.pos_ = pbc.apply(particle.xform_.pos_);
        particle.time = time;
    }
}

int main(int argc, char *argv[]){


    Configuration config;
    xml_load_config(argv[1], config);

    Simulation sim(config);

    PeriodicCallback output(0.001);
    output.setNextFunction([](double time){
        return time + 0.1;
    });
    output.setCallback([&sim](double time){
        static int nFiles = 0;
        Configuration config = sim.get_configuration();
        double kinetic = 0.0;
        for(auto& particle: config.particles_){
            update_particle(particle, time, config.pbc_);
            kinetic += particle.vel.length2();
        }
        char buff[64];
        sprintf(buff, "Data/pid%u.%06u.xml", getpid(), ++nFiles);
        xml_save_config(buff, config);
        clam::Vec3d box_size = config.pbc_.getSize();
        double volume = box_size[0] * box_size[1] * box_size[2];
        double kT = kinetic / (3.0 * config.particles_.size());
        double pressure = (config.particles_.size() - sim.get_stress() / (3.0 * time * kT)) / volume;
        printf("%f\t%f\n", pressure, kT);
    });

    sim.run(1000.0, output);

    return 0;
}
