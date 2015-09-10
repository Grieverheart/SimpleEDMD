#include <cstdio>
#include <cstring>
#include <unistd.h>
#include "shape/variant.h"
#include "Simulation.h"

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
                if(u < 3)	coords[u] = atof(t1);
                else radius = atof(t1);
                u++;
            }
            Particle part;
            part.pos = pbc.apply(coords);
            part.rot = clam::Quatd(1.0, 0.0, 0.0, 0.0);
            part.radius = radius;
            particles.push_back(part);
		}
        else if(i == 2) pbc.setSize(atof(line));
		i++;
	}
	free(t1);
	fclose(fp);

    shapes.push_back(new shape::Variant(shape::Sphere()));
}

void saveConfig(const char* filename, double time, const Simulation& sim){
    FILE *fp = fopen(filename, "w");
    int nSpheres = sim.getNumParticles();
    const CubicPBC& pbc = sim.getPBC();
    fprintf(fp, "%d\n%f\n", nSpheres, pbc.getSize());

    const std::vector<Particle>& particles = sim.getParticles();
    for(int i = 0; i < nSpheres; ++i){
        clam::Vec3d pos(particles[i].pos);
        pos += particles[i].vel * (time - particles[i].time) - sim.getSystemVelocity() * time;
        pos  = pbc.apply(pos);
        fprintf(fp, "%f\t%f\t%f\t", pos[0], pos[1], pos[2]);
        fprintf(fp, "%f\n", particles[i].radius);
    }
    fprintf(fp, "%f\n", time);
    fclose(fp);
}

int main(int argc, char *argv[]){


    std::vector<Particle> particles;
    std::vector<shape::Variant*> shapes;
    CubicPBC pbc;
    readConfig(argv[1], pbc, particles, shapes);

    Simulation sim(pbc, std::move(particles), std::move(shapes));

    sim.init();

    PeriodicCallback output(100.0 + 0.001);
    output.setNextFunction([](double time){
        return 100.0 + 1.175 * (time - 100.0);
    });
    output.setCallback([&sim](double time){
        static int nFiles = 0;
        char buff[64];
        sprintf(buff, "Data/pid%u.%06u.dat", getpid(), nFiles++);
        saveConfig(buff, time, sim);
    });

    sim.run(500.0, output);

    return 0;
}
