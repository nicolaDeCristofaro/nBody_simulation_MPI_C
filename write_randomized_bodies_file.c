#include <stdio.h>
#include <stdlib.h>

typedef struct {
    float mass;  //massa della particella
    float x, y, z; //posizione della particella identificata dalle coordinate x,y,z
    float vx, vy, vz; //velocità della particella nelle rispettive direzioni x,y,z ??
} Particle;

/*Function that initializes the collection of particles*/
void randomizeParticles(Particle *particles, int n) {
    for (int i = 0; i < n; i++) {
        particles[i].mass = 1.0;  // valore arbitrario scelto per la massa delle particelle -> 1.0
        particles[i].x = 2.0f * (rand() / (float)RAND_MAX) - 1.0f;  //random number between -1 and 1
		particles[i].y = 2.0f * (rand() / (float)RAND_MAX) - 1.0f;
		particles[i].z = 2.0f * (rand() / (float)RAND_MAX) - 1.0f;
        //particles[i].vx = 2.0f * (rand() / (float)RAND_MAX) - 1.0f;
		//particles[i].vy = 2.0f * (rand() / (float)RAND_MAX) - 1.0f;
		//particles[i].vz = 2.0f * (rand() / (float)RAND_MAX) - 1.0f;
  }
}


int main(int argc, char* argv[]){


    int numParticles = 1000;  //default number of particles if no parameters is given from command-line
    //Prendo da linea di comando il numero di elementi dell'array
    if(argc > 1){
        //il programma ha 2 argomenti: 1 standard(nome del programma) e 2 quello passato da linea di comando
        //questa è la configurazione che ci serve
        numParticles = atoi(argv[1]);
    }

    int bytes = numParticles*sizeof(Particle);
    Particle *particles = (Particle*)malloc(bytes);

    srand (0);
    randomizeParticles(particles, numParticles);

    FILE * file= fopen("output.txt", "wb");
    if (file != NULL) {
        fwrite(particles , sizeof(Particle)*numParticles, 1, file);
        fclose(file);
    }

    Particle *particles2 = (Particle*)malloc(bytes);
    file= fopen("output.txt", "rb");
    int i=0;
    if (file != NULL) {
        fread(particles2, sizeof(Particle)*numParticles, 1, file);
        fclose(file);
    }

    for(int i=0; i< numParticles; i++){
        printf("[%d].x = %.2f\t", i, particles2[i].x);
        printf("[%d].y = %.2f\t", i, particles2[i].y);
        printf("[%d].z = %.2f\t", i, particles2[i].z);
        printf("[%d].vx = %.2f\t", i, particles2[i].vx);
        printf("[%d].vy = %.2f\t", i, particles2[i].vy);
        printf("[%d].vz = %.2f\t", i, particles2[i].vz);
        printf("\n");
    }


}