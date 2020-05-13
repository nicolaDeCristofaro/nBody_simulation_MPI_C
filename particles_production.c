#include <stdio.h>
#include <stdlib.h>
#include <math.h>

typedef struct{
    float mass;
    float x, y, z;
    float vx, vy, vz;
} Particle;

/*Inizializza con valori random lo stato delle particelle*/
void randomizeParticles(Particle *particles, int n) {
    for (int i = 0; i < n; i++) {
        particles[i].mass = 2.0;  // valore arbitrario scelto per la massa delle particelle -> 2.0
        particles[i].x = 2.0f * (rand() / (float)RAND_MAX) - 1.0f;  //numero random compreso tra -1 and 1
		particles[i].y = 2.0f * (rand() / (float)RAND_MAX) - 1.0f;
		particles[i].z = 2.0f * (rand() / (float)RAND_MAX) - 1.0f;
        particles[i].vx = 2.0f * (rand() / (float)RAND_MAX) - 1.0f;
	    particles[i].vy = 2.0f * (rand() / (float)RAND_MAX) - 1.0f;
		particles[i].vz = 2.0f * (rand() / (float)RAND_MAX) - 1.0f;
  }
}

int main(int argc, char* argv[]){

    int num_particles = 1000;  //Numero delle particelle di DEFAULT se nessun parametro è fornito sulla command-line
    if(argc > 1){
        // E' stato fornito il parametro da linea di comando che indica il numero di particelle
        num_particles = atoi(argv[1]);
    }

    Particle *particles = NULL;
    particles = (Particle*) malloc(num_particles * sizeof(Particle));

    randomizeParticles(particles,num_particles);

    FILE * file= fopen("particles.dat", "w");
    if (file != NULL) {
        fwrite(particles , sizeof(Particle) * num_particles, 1, file);
        fclose(file);
    }

    free(particles);

    /* TEST: decommenta per scrivere su stdout il valore iniziale delle particelle
    printf("INPUT\n");
    for (int i = 0  ; i < num_particles; i++) { 
        printf("[%d].x = %f\t", i, particles[i].x);
        printf("[%d].y = %f\t", i, particles[i].y);
        printf("[%d].z = %f\t", i, particles[i].z);
        printf("[%d].vx = %f\t", i, particles[i].vx);
        printf("[%d].vy = %f\t", i, particles[i].vy);
        printf("[%d].vz = %f\t", i, particles[i].vz);
        printf("\n");
    }*/
}