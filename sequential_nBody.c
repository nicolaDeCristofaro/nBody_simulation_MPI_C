#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define SOFTENING 1e-9f

/*Implementing the simulation of the n-body problem
  Sequential_version*/

typedef struct
{
  float mass;
  float x, y, z;
  float vx, vy, vz;
} Particle;

/*Inizializza con valori random lo stato delle particelle*/
void randomizeParticles(Particle *particles, int n){
  for (int i = 0; i < n; i++){
    particles[i].mass = 2.0;                                   // valore arbitrario scelto per la massa delle particelle -> 2.0
    particles[i].x = 2.0f * (rand() / (float)RAND_MAX) - 1.0f; //numero random compreso tra -1 and 1
    particles[i].y = 2.0f * (rand() / (float)RAND_MAX) - 1.0f;
    particles[i].z = 2.0f * (rand() / (float)RAND_MAX) - 1.0f;
    particles[i].vx = 2.0f * (rand() / (float)RAND_MAX) - 1.0f;
    particles[i].vy = 2.0f * (rand() / (float)RAND_MAX) - 1.0f;
    particles[i].vz = 2.0f * (rand() / (float)RAND_MAX) - 1.0f;
  }
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

int main(const int argc, const char **argv){

  int nBodies = 10; //number of bodies if no paramters is given from command-line
  if (argc > 1)
    nBodies = atoi(argv[1]);

  const float dt = 0.01f; // time step
  const int nIters = 10;  // simulation iterations

  int bytes = nBodies * sizeof(Particle);
  float *buf = (float *)malloc(bytes);
  Particle *p = (Particle *)buf;

  srand(0);
  randomizeParticles(p, nBodies);

  /*INPUT TEST
  printf("INPUT\n");
  for(int i=0; i< nBodies; i++){
    printf("[%d].x = %f  ", i, p[i].x);
    printf("[%d].y = %f  ", i, p[i].y);
    printf("[%d].z = %f  ", i, p[i].z);
    printf("[%d].vx = %f  ", i, p[i].vx);
    printf("[%d].vy = %f  ", i, p[i].vy);
    printf("[%d].vz = %f  ", i, p[i].vz);
    printf("\n");
  }*/

  clock_t start;
  clock_t end;
  double totalTime = 0.0;

  for (int iter = 1; iter <= nIters; iter++){
    start = clock();

    bodyForce(p, dt, nBodies); // compute interbody forces

    for (int i = 0; i < nBodies; i++){ // integrate position
      p[i].x += p[i].vx * dt;
      p[i].y += p[i].vy * dt;
      p[i].z += p[i].vz * dt;
    }

    end = clock() - start;

    printf("Iteration %d: %f seconds\n", iter, (double)end / CLOCKS_PER_SEC);
    totalTime += (double)end / CLOCKS_PER_SEC;
  }

  double avgTime = totalTime / (double)(nIters - 1);
  printf("Total time: %f seconds\n", totalTime);
  printf("Avg time: %f seconds\n", avgTime);

  /*OUTPUT TEST
  printf("OUTPUT\n");
  for (int i = 0; i < nBodies; i++)
  {
    printf("[%d].x = %f  ", i, p[i].x);
    printf("[%d].y = %f  ", i, p[i].y);
    printf("[%d].z = %f  ", i, p[i].z);
    printf("[%d].vx = %f  ", i, p[i].vx);
    printf("[%d].vy = %f  ", i, p[i].vy);
    printf("[%d].vz = %f  ", i, p[i].vz);
    printf("\n");
  }*/

  /*Scrivo l'output su file per poi poterne valutare la correttenza confrontando con l'output parallelo*/
  FILE *file = fopen("sequential_output.txt", "wb");
  if (file != NULL){
    fwrite(p, sizeof(Particle) * nBodies, 1, file);
    fclose(file);
  }

  free(buf);
}
