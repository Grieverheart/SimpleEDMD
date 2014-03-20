#include <cstdio>
#include <cstring>
#include "include/Simulation.h"

void readConfig(const char* filename, CubicPBC& pbc, std::vector<Particle>& particles){
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
            Vec3d coords;
            double radius = 0.0;
            for(t1 = strtok(line,delim); t1 != NULL; t1 = strtok(NULL, delim)){
                if(u < 3)	coords[u] = atof(t1);
                else radius = atof(t1);
                u++;
            }
            Particle part;
            part.pos = pbc.apply(coords);
            part.radius = radius;
            particles.emplace_back(part);
		}
        else if(i == 2) pbc.setSize(atof(line));
		i++;
	}
	free(t1);
	fclose(fp);
}

int main(int argc, char *argv[]){


    std::vector<Particle> particles;
    CubicPBC pbc;
    readConfig(argv[1], pbc, particles);

    Simulation sim(pbc, std::move(particles));

    sim.init();

    PeriodicCondition end(100.0);
    PeriodicCondition output(200.0);
    output.setNextFunction([](double time){
        return time + 100.0;
    });
    sim.run(output, end);

    return 0;
}
