#include <cstdio>
#include <cstring>
#include <unistd.h>
#include "shape/variant.h"
#include "Simulation.h"
#include "obj_loader.h"

#if 1
void readConfig(const char* filename, CubicPBC& pbc, std::vector<Particle>& particles, std::vector<shape::Variant*>& shapes){
    char line[128];
    const char* delim = "\t";
    char *t1 = NULL;
    int i = 1,u = 0;

    FILE *fp;
    fp = fopen(filename,"r");
    if(!fp){
        printf("Couldn't open file %s\n",filename);
    }

    while(fgets(line,sizeof(line),fp) != NULL){
        u=0;
        if(line[strlen(line)-1] == '\n') line[strlen(line)-1] = '\0'; // Remove the newline from the end of the string
        if(i > 2){
            Particle part;
            for(t1 = strtok(line,delim); t1 != NULL; t1 = strtok(NULL, delim)){
                if(u < 3) part.pos[u] = atof(t1);
                else part.rot[u - 3] = atof(t1);
                u++;
            }
            part.pos = pbc.apply(part.pos);
            part.size = 1.0;
            part.shape_id = 0;
            particles.push_back(part);
        }
        else if(i == 2) pbc.setSize(atof(line));
        ++i;
    }
    free(t1);
    fclose(fp);

    std::vector<clam::Vec3d> vertices;
    std::vector<std::vector<unsigned int>> faces;
    if(!load_obj("obj/Rhombic123.obj", vertices, faces)) printf("Couldn't load Rhombic123.obj!\n");
    shapes.push_back(new shape::Variant(shape::Polyhedron(vertices, faces)));
}

#else
void readConfig(const char* filename, CubicPBC& pbc, std::vector<Particle>& particles, std::vector<shape::Variant*>& shapes){
    char line[128];
    const char* delim = "\t";
    char *t1 = NULL;
    int i = 1,u = 0;

    FILE *fp;
    fp = fopen(filename,"r");
    if(!fp){
        printf("Couldn't open file %s\n",filename);
    }

    while(fgets(line,sizeof(line),fp) != NULL){
        u=0;
        if(line[strlen(line)-1] == '\n') line[strlen(line)-1] = '\0'; // Remove the niewline from the end of the string
        if(i > 2){
            clam::Vec3d coords;
            double radius = 0.0;
            for(t1 = strtok(line,delim); t1 != NULL; t1 = strtok(NULL, delim)){
                if(u < 3) coords[u] = atof(t1);
                else radius = atof(t1);
                u++;
            }
            Particle part;
            part.pos = pbc.apply(coords);
            part.rot = clam::Quatd(1.0, 0.0, 0.0, 0.0);
            part.size = radius;
            part.shape_id = 0;
            particles.push_back(part);
        }
        else if(i == 2) pbc.setSize(atof(line));
        i++;
    }
    free(t1);
    fclose(fp);

    shapes.push_back(new shape::Variant(shape::Sphere()));
}
#endif

void saveConfig(const char* filename, double time, const Simulation& sim){
    FILE *fp = fopen(filename, "w");
    size_t nSpheres = sim.get_particles().size();
    const CubicPBC& pbc = sim.get_pbc();
    fprintf(fp, "%lu\n", nSpheres);
    fprintf(fp, "%f\t0.0\t0.0\t0.0\t%f\t0.0\t0.0\t0.0\t%f\n", pbc.getSize(), pbc.getSize(), pbc.getSize());

    const std::vector<Particle>& particles = sim.get_particles();
    for(size_t i = 0; i < nSpheres; ++i){
        clam::Vec3d pos(particles[i].pos);
        pos += particles[i].vel * (time - particles[i].time) - sim.get_system_velocity() * time;
        pos  = pbc.apply(pos);
        clam::Quatd rot = particles[i].rot;
        {
            double ang_vel_abs = particles[i].ang_vel.length();
            if(ang_vel_abs != 0.0){
                rot = clam::fromAxisAngle(
                    ang_vel_abs * (time - particles[i].time),
                    particles[i].ang_vel / ang_vel_abs
                ) * rot;
            }
        }
        fprintf(fp, "%f\t%f\t%f\t", pos[0], pos[1], pos[2]);
        clam::Vec3d axis;
        double angle;
        rot.toAxisAngle(angle, axis);
        fprintf(fp, "%f\t%f\t%f\t%f\n", angle * 180.0 / M_PI, axis[0], axis[1], axis[2]);
    }
    //fprintf(fp, "%f\n", time);
    fclose(fp);
}

int main(int argc, char *argv[]){


    std::vector<Particle> particles;
    std::vector<shape::Variant*> shapes;
    CubicPBC pbc;
    readConfig(argv[1], pbc, particles, shapes);

    Simulation sim(Configuration(pbc, particles, shapes));

    char buff[64];
    sprintf(buff, "Data/pid%u.%06u.dat", getpid(), 0);
    saveConfig(buff, 0.0, sim);

    PeriodicCallback output(0.001);
    output.setNextFunction([](double time){
        return time + 0.1;
    });
    output.setCallback([&sim](double time){
        static int nFiles = 0;
        char buff[64];
        sprintf(buff, "Data/pid%u.%06u.dat", getpid(), ++nFiles);
        saveConfig(buff, time, sim);
    });

    sim.run(1.0, output);

    return 0;
}
