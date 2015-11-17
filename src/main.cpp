#include <cstdio>
#include <cstring>
#include <unistd.h>
#include "shape/variant.h"
#include "Simulation.h"
#include "obj_loader.h"
#include "io/config_xml.h"

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
        for(auto& particle: config.particles_){
            update_particle(particle, time, config.pbc_);
        }
        char buff[64];
        sprintf(buff, "Data/pid%u.%06u.xml", getpid(), ++nFiles);
        xml_save_config(buff, config);
    });

    sim.run(1.0, output);

    return 0;
}
