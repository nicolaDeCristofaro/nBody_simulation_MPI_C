#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define SOFTENING 1e-9f

/*Implementing the simulation of the n-body problem
  Sequential_version, taken from the example given 
  at link https://github.com/harrism/mini-nbody/blob/master/nbody.c
  and adapted*/

typedef struct
{
  float mass;
  float x, y, z;
  float vx, vy, vz;
} Particle;

//Function definition
void bodyForce(Particle *p, float dt, int n);


int main(const int argc, const char **argv){

  int nBodies = 10; //number of bodies if no paramters is given from command-line
  if (argc > 1)
    nBodies = atoi(argv[1]);

  const float dt = 0.01f; // time step
  const int nIters = 10;  // simulation iterations

  Particle *particles = NULL;
  particles = (Particle*) malloc(nBodies * sizeof(Particle));

  FILE *fileRead = fopen("particles.txt", "r");
  if (fileRead == NULL){
      /* Impossibile aprire il file */
      printf("\nImpossibile aprire il file.\n");
      exit(EXIT_FAILURE);
  }

  fread(particles, sizeof(Particle) * nBodies, 1, fileRead);
  fclose(fileRead);

  clock_t start;
  clock_t end;
  double totalTime = 0.0;

  for (int iter = 1; iter <= nIters; iter++){
    start = clock();

    bodyForce(particles, dt, nBodies); // compute interbody forces

    for (int i = 0; i < nBodies; i++){ // integrate position
      particles[i].x += particles[i].vx * dt;
      particles[i].y += particles[i].vy * dt;
      particles[i].z += particles[i].vz * dt;
    }

    end = clock() - start;

    //printf("Iteration %d: %f seconds\n", iter, (double)end / CLOCKS_PER_SEC);
    totalTime += (double)end / CLOCKS_PER_SEC;
  }

  double avgTime = totalTime / (double)(nIters - 1);
  printf("\nComputation completed\n");
  printf("Avg iteration time: %f seconds\n", avgTime);
  printf("Total time: %f\n", totalTime);

  free(particles);
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