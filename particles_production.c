#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <errno.h>
#include <limits.h>

typedef struct{
    float mass;
    float x, y, z;
    float vx, vy, vz;
} Particle;

/*Definizione di funzioni*/
int convertStringToInt(char *str);
void randomizeParticles(Particle *particles, int n);

int main(int argc, char* argv[]){

    int num_particles = 1000;  //Numero delle particelle di DEFAULT se nessun parametro è fornito sulla command-line
    if(argc > 1){
        // E' stato fornito il parametro da linea di comando che indica il numero di particelle
        //Converti stringa a intero
        num_particles = convertStringToInt(argv[1]);
    }

    Particle *particles = NULL;
    particles = (Particle*) malloc(num_particles * sizeof(Particle));

    srand(0);
    randomizeParticles(particles,num_particles);

    FILE * file= fopen("particles.txt", "w");
    if (file != NULL) {
        fwrite(particles , sizeof(Particle) * num_particles, 1, file);
        fclose(file);
        printf("%d particelle sono state create con valori random e scritte su file\n", num_particles);
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