#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <limits.h>
#include <errno.h>

#define SOFTENING 1e-9f

/*Implementing the simulation of the n-body problem
  Sequential_version, taken from the example given 
  at link https://github.com/harrism/mini-nbody/blob/master/nbody.c
  and adapted*/

typedef struct{
    float mass;
    float x, y, z;
    float vx, vy, vz;
} Particle;

/* Functions definition */
int convertStringToInt(char *str);
void bodyForce(Particle *p, float dt, int n);

int main(int argc, char* argv[]){

    int nBodies = 10; //number of bodies if no paramters is given from command-line
    if (argc > 1) nBodies = convertStringToInt(argv[1]);

    const float dt = 0.01f; // time step
    const int nIters = 10;  // simulation iterations

    clock_t startIter, endIter;
    clock_t startTotal = clock(), endTotal;
    double totalTime = 0.0;


    Particle *particles = NULL;
    particles = (Particle *)malloc(nBodies * sizeof(Particle));

    FILE *fileRead = fopen("particles.txt", "r");
    if (fileRead == NULL){
        /* Impossibile aprire il file */
        printf("\nImpossibile aprire il file.\n");
        exit(EXIT_FAILURE);
    }

    int particlesRead = fread(particles, sizeof(Particle) * nBodies, 1, fileRead);
    if( particlesRead == 0){
        /*il numero di particelle da leggere è maggiore del numero di particelle nel file*/
        printf("ERROR: Il numero di particelle da leggere è maggiore del numero di particelle nel file\n");
        exit(EXIT_FAILURE);
    }
    fclose(fileRead);

    /* TEST: decommenta per scrivere su stdout lo stato iniziale delle particelle dopo la lettura da file
    printf("INPUT\n");
    for(int i=0; i< nBodies; i++){
        printf("[%d].x = %f\t", i, particles[i].x);
        printf("[%d].y = %f\t", i, particles[i].y);
        printf("[%d].z = %f\t", i, particles[i].z);
        printf("[%d].vx = %f\t", i, particles[i].vx);
        printf("[%d].vy = %f\t", i, particles[i].vy);
        printf("[%d].vz = %f\t", i, particles[i].vz);
        printf("\n");
    }*/

    for (int iter = 1; iter <= nIters; iter++){
        startIter = clock();

        bodyForce(particles, dt, nBodies); // compute interbody forces

        for (int i = 0; i < nBodies; i++){ // integrate position
            particles[i].x += particles[i].vx * dt;
            particles[i].y += particles[i].vy * dt;
            particles[i].z += particles[i].vz * dt;
        }

        endIter = clock() - startIter;
        printf("Iterazione %d di %d completata in %f seconds\n", iter, nIters, (double)endIter / CLOCKS_PER_SEC);
    }

    endTotal = clock();
    totalTime = (double)(endTotal - startTotal) / CLOCKS_PER_SEC;
    double avgTime = totalTime / (double)(nIters);
    printf("\nAvg iteration time: %f seconds\n", avgTime);
    printf("Total time: %f seconds\n", totalTime);

    /* TEST: decommenta per scrivere su stdout lo stato finale delle particelle dopo la computazione
    printf("OUTPUT\n");
    for (int i = 0; i < nBodies; i++){
        printf("[%d].x = %f  ", i, particles[i].x);
        printf("[%d].y = %f  ", i, particles[i].y);
        printf("[%d].z = %f  ", i, particles[i].z);
        printf("[%d].vx = %f  ", i, particles[i].vx);
        printf("[%d].vy = %f  ", i, particles[i].vy);
        printf("[%d].vz = %f  ", i, particles[i].vz);
        printf("\n");
    }*/

    /*Scrivo l'output su file per poi poterne valutare la correttenza confrontando con l'output parallelo*/
    FILE *fileWrite = fopen("sequential_output.txt", "w");
    if (fileWrite != NULL){
        fwrite(particles, sizeof(Particle) * nBodies, 1, fileWrite);
        fclose(fileWrite);
    }

    free(particles);

}

/*Conversione da stringa a intero*/
int convertStringToInt(char *str){
    char *endptr;
    long val;  
    errno = 0;  //Per distinguere successo/fallimento dopo la chiamata

    val = strtol(str, &endptr, 10);

    /* Check possibili errori */
    if ((errno == ERANGE && (val == LONG_MAX || val == LONG_MIN)) || (errno != 0 && val == 0)) {
        perror("strtol");
        exit(EXIT_FAILURE);
    }

    if (endptr == str) {
        fprintf(stderr, "No digits were found\n");
        exit(EXIT_FAILURE);
    }

    /* Se siamo qui, strtol() ha convertito un numero correttamente */
    return (int)val;
}

/*Function that make computation*/
void bodyForce(Particle *p, float dt, int n){

    for (int i = 0; i < n; i++){
        float Fx = 0.0f;
        float Fy = 0.0f;
        float Fz = 0.0f;

        for (int j = 0; j < n; j++){
            float dx = p[j].x - p[i].x;
            float dy = p[j].y - p[i].y;
            float dz = p[j].z - p[i].z;
            float distSqr = dx * dx + dy * dy + dz * dz + SOFTENING;
            float invDist = 1.0f / sqrtf(distSqr);
            float invDist3 = invDist * invDist * invDist;

            Fx += dx * invDist3;
            Fy += dy * invDist3;
            Fz += dz * invDist3;
        }

        p[i].vx += dt * Fx;
        p[i].vy += dt * Fy;
        p[i].vz += dt * Fz;
    }
}

  