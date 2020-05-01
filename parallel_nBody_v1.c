#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#define MASTER 0            // Rank of the Master processor  
#define I 10               // Number of iterations
#define SOFTENING 1e-9f     // Infinitely large value     


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
        particles[i].vx = 2.0f * (rand() / (float)RAND_MAX) - 1.0f;
		particles[i].vy = 2.0f * (rand() / (float)RAND_MAX) - 1.0f;
		particles[i].vz = 2.0f * (rand() / (float)RAND_MAX) - 1.0f;
  }
}

/*Funzione che esegue computazione: calcola le nuove posizioni e velocità di ogni particella 
che dipendono dalla forza esercitata da tutte le altre particelle (complessità -> n^2) */
void bodyForce(Particle *particles, float dt, int n) {
    for (int i = 0; i < n; i++) { 
        float Fx = 0.0f; float Fy = 0.0f; float Fz = 0.0f;

        for (int j = 0; j < n; j++) {
            float dx = particles[j].x - particles[i].x;
            float dy = particles[j].y - particles[i].y;
            float dz = particles[j].z - particles[i].z;
            float distSqr = dx*dx + dy*dy + dz*dz + SOFTENING;
            float invDist = 1.0f / sqrtf(distSqr);
            float invDist3 = invDist * invDist * invDist;

            Fx += dx * invDist3; Fy += dy * invDist3; Fz += dz * invDist3;
        }

    particles[i].vx += dt*Fx; 
    particles[i].vy += dt*Fy; 
    particles[i].vz += dt*Fz;
    }
}

int main(int argc, char* argv[]){

    int  numtasks, myrank;
    double start, end, iterStart, iterEnd;

    /***** MPI Initializations *****/
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
    MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
    MPI_Barrier(MPI_COMM_WORLD); /* tutti i processi sono inizializzati */
    start = MPI_Wtime();


    int numParticles = 1000;  //default number of particles if no parameters is given from command-line
    //Prendo da linea di comando il numero di elementi dell'array
    if(argc > 1){
        //il programma ha 2 argomenti: 1 standard(nome del programma) e 2 quello passato da linea di comando
        //questa è la configurazione che ci serve
        numParticles = atoi(argv[1]);
    }

    const float dt = 0.01f; // time step

    int bytes = numParticles*sizeof(Particle);
    Particle *particles = (Particle*)malloc(bytes);

    randomizeParticles(particles, numParticles);

    /*printf("Print the array of bodies initialized\n");
    for(int i=0; i< numParticles; i++){
        printf("[%d].x = %.2f\n", i, particles[i].x);
        printf("[%d].y = %.2f\n", i, particles[i].y);
        printf("[%d].z = %.2f\n", i, particles[i].z);
        printf("[%d].vx = %.2f\n", i, particles[i].vx);
        printf("[%d].vy = %.2f\n", i, particles[i].vy);
        printf("[%d].vz = %.2f\n", i, particles[i].vz);
    }*/

    for (int iter = 0; iter < I; iter++) {

        // Synchronize before starting timing
        MPI_Barrier(MPI_COMM_WORLD);
        iterStart = MPI_Wtime();

        bodyForce(particles, dt, numParticles); // compute interbody forces

        for (int i = 0 ; i < numParticles; i++) { // integrate position
            particles[i].x += particles[i].vx*dt;
            particles[i].y += particles[i].vy*dt;
            particles[i].z += particles[i].vz*dt;
        }

        MPI_Barrier(MPI_COMM_WORLD);
        iterEnd = MPI_Wtime();
        printf("Iteration %d: %f seconds\n", iter, iterEnd-iterStart);

    }

    MPI_Barrier(MPI_COMM_WORLD); /* tutti i processi hanno terminato */
    end = MPI_Wtime();
    MPI_Finalize();

    double totalTime = end-start;
    double avgTime = totalTime / (double)(I); 
    printf("Avg time: %f seconds\n", avgTime);
    printf("Total time: %f\n", totalTime);

    free(particles);

}
