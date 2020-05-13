#include <stdio.h>
#include <stdlib.h>

typedef struct{
    float mass;
    float x, y, z;
    float vx, vy, vz;
} Particle;

//confronta due float e ritorna 1 se sono uguali 0 altrimenti
int compare_float(float f1, float f2){
    float precision = 0.000001;
    if (((f1 - precision) < f2) &&
        ((f1 + precision) > f2))
        return 1;
    else
        return 0;
}

/*Programma per valutare la correttezza dell'implementazione parallela:
Confronta l'output della versione sequenziale con l'output della soluzione parallela*/
int main(int argc, char *argv[]){

    FILE *seq_file = fopen("sequential_output.txt", "rb");
    FILE *par_file = fopen("parallel_output.txt", "rb");

    int num_particles = 5;

    /* fopen() ritorna NULL se non riesce ad aprire un file nella modalità indicata. */
    if (seq_file == NULL || par_file == NULL){
        /* Impossibile aprire uno dei file quindi esco */
        printf("\nImpossibile aprire uno dei files.\n");
        printf("Controlla se i files esistono e se hai i privilegi per leggerli.\n");
        exit(EXIT_FAILURE);
    }

    Particle *seq_particles = (Particle *)malloc(num_particles * sizeof(Particle));
    Particle *par_particles = (Particle *)malloc(num_particles * sizeof(Particle));

    //Leggo lo stato finale delle particelle dalla computazione sequenziale e da quella parallela
    fread(seq_particles, sizeof(Particle) * num_particles, 1, seq_file);
    fread(par_particles, sizeof(Particle) * num_particles, 1, par_file);

    /*Verifico se lo stato finale delle particelle dopo la computazione sequenziale 
    è uguale a quello delle particelle dopo la computazione parallela*/
    for (int i = 0; i < num_particles; i++){
        if (!compare_float(seq_particles[i].x, par_particles[i].x))
            printf("Difference in member [%d].x \n", i);
        if (!compare_float(seq_particles[i].y, par_particles[i].y))
            printf("Difference in member [%d].y \n", i);
        if (!compare_float(seq_particles[i].z, par_particles[i].z))
            printf("Difference in member [%d].z \n", i);
        if (!compare_float(seq_particles[i].vx, par_particles[i].vx))
            printf("Difference in member [%d].vx \n", i);
        if (!compare_float(seq_particles[i].vy, par_particles[i].vy))
            printf("Difference in member [%d].vy \n", i);
        if (!compare_float(seq_particles[i].vz, par_particles[i].vz))
            printf("Difference in member [%d].vz \n", i);
    }

    /* Infine chiudo i files per rilasciare le risorse */
    fclose(seq_file);
    fclose(par_file);
    free(seq_particles);
    free(par_particles);

    return 0;
}