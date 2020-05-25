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

int main(const int argc, const char **argv)
{

  int nBodies = 10; //number of bodies if no paramters is given from command-line
  if (argc > 1)
    nBodies = atoi(argv[1]);

  const float dt = 0.01f; // time step
  const int nIters = 10;  // simulation iterations

  clock_t startIter, endIter;
  clock_t startTotal = clock(), endTotal;
  double totalTime = 0.0;


  Particle *particles = NULL;
  particles = (Particle *)malloc(nBodies * sizeof(Particle));

  /*INPUT TEST
  printf("INPUT\n");
  for(int i=0; i< nBodies; i++){
    printf("[%d].x = %f  ", i, particles[i].x);
    printf("[%d].y = %f  ", i, particles[i].y);
    printf("[%d].z = %f  ", i, particles[i].z);
    printf("[%d].vx = %f  ", i, particles[i].vx);
    printf("[%d].vy = %f  ", i, particles[i].vy);
    printf("[%d].vz = %f  ", i, particles[i].vz);
    printf("\n");
  }*/

  FILE *fileRead = fopen("particles.txt", "r");
  if (fileRead == NULL)
  {
    /* Impossibile aprire il file */
    printf("\nImpossibile aprire il file.\n");
    exit(EXIT_FAILURE);
  }

  fread(particles, sizeof(Particle) * nBodies, 1, fileRead);
  fclose(fileRead);


  for (int iter = 1; iter <= nIters; iter++)
  {
    startIter = clock();

    bodyForce(particles, dt, nBodies); // compute interbody forces

    for (int i = 0; i < nBodies; i++)
    { // integrate position
      particles[i].x += particles[i].vx * dt;
      particles[i].y += particles[i].vy * dt;
      particles[i].z += particles[i].vz * dt;
    }

    endIter = clock() - startIter;
    printf("Iterazione %d di %d completata in %f seconds\n", iter, nIters, (double)endIter / CLOCKS_PER_SEC);
  }

  endTotal = clock();
  totalTime = (double)(endTotal - startTotal) / CLOCKS_PER_SEC;
  double avgTime = totalTime / (double)(nIters - 1);
  printf("\nAvg iteration time: %f seconds\n", avgTime);
  printf("Total time: %f seconds\n", totalTime);

  /*OUTPUT TEST
  printf("OUTPUT\n");
  for (int i = 0; i < nBodies; i++)
  {
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
  if (fileWrite != NULL)
  {
    fwrite(particles, sizeof(Particle) * nBodies, 1, fileWrite);
    fclose(fileWrite);
  }

  free(particles);
}

/*Function that make computation*/
void bodyForce(Particle *p, float dt, int n)
{
  for (int i = 0; i < n; i++)
  {
    float Fx = 0.0f;
    float Fy = 0.0f;
    float Fz = 0.0f;

    for (int j = 0; j < n; j++)
    {
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
