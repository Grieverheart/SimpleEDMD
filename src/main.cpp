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
    const double output_delta = 0.1;
    double output_start_time = 0.01;

    const auto& directory = argv[3];

    if(strcmp(argv[1], "-r") == 0){
        FILE* fp = fopen(argv[2], "rb");
        fseek(fp, 0l, SEEK_END);
        size_t file_size = ftell(fp);
        fseek(fp, 0l, SEEK_SET);
        char* data = reinterpret_cast<char*>(malloc(file_size));
        fread(data, file_size, 1, fp);
        fclose(fp);

        Archive ar(data, file_size);
        sim = new Simulation();
        deserialize(ar, sim);

        free(data);

        output_start_time = sim->time() + output_delta;
    }
    else{
        Configuration config;
        xml_load_config(argv[1], config);

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

    double pf = packing_fraction(sim->configuration());

    PeriodicCallback output(output_start_time);
    output.setNextFunction([output_delta](double time){
        return time + output_delta;
    });

    FILE* pressure_fp;
    {
        char fp_buff[64];
        sprintf(fp_buff, "%s/pressure.pf%.3f.pid%u.dat", directory, pf, getpid());
        pressure_fp = fopen(fp_buff, "w");
    }

    output.setCallback([sim, pf, pressure_fp, directory](double time){
        printf("%f\n", time);
        static int nFiles = 0;
        Configuration config = sim->configuration();

        for(auto& particle: config.particles_){
            update_particle(particle, time, config.pbc_);
        }

        char buff[64];
        sprintf(buff, "%s/pid%u.pf%.3f.step%06u.xml", directory, getpid(), pf, ++nFiles);
        xml_save_config(buff, config);

        fprintf(pressure_fp, "%e\t%e\n", sim->time(), sim->average_pressure());
        fflush(pressure_fp);

        sim->reset_statistics();

//#ifndef NDEBUG
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
        //assert(overlaps == false);
        //if(overlaps) exit(0);
//#endif
        Archive ar;
        serialize(ar, *sim);
        sprintf(buff, "%s/archive.pf%.3f.pid%u.bin", directory, pf, getpid());
        FILE* fp = fopen(buff, "wb");
        fwrite(ar.data(), 1, ar.size(), fp);
        fclose(fp);
    });

    sim->run(1.0, output);

    fclose(pressure_fp);

    return 0;
}
