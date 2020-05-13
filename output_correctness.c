#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>

typedef struct{
    float mass;
    float x, y, z;
    float vx, vy, vz;
} Particle;


//confronta due float e ritorna 1 se sono uguali 0 altrimenti
int compare_float(float f1, float f2){

    float precision = 0.01f;
    //if( (fabs(f1 - f2)) <= precision) return 1;
    //else return 0;
    if (((f1 - precision) < f2) &&
        ((f1 + precision) > f2))
        return 1;
    else
        return 0;
}

/*Programma per valutare la correttezza dell'implementazione parallela:
Confronta l'output della versione sequenziale con l'output della soluzione parallela*/
int main(int argc, char *argv[]){

    FILE *seq_file = fopen("sequential_output.txt", "r");
    FILE *par_file = fopen("parallel_output.txt", "r");

    int num_particles = 1000;  //Numero delle particelle di DEFAULT se nessun parametro è fornito sulla command-line
    if(argc > 1){
        // E' stato fornito il parametro da linea di comando che indica il numero di particelle
        num_particles = atoi(argv[1]);
    }

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
    int flag = 1; //flag rimarrà uguale a 1 se i due autput sono uguali altrimenti sarà 0
    for (int i = 0; i < num_particles; i++){
        if (!compare_float(seq_particles[i].x, par_particles[i].x)){
            printf("Difference in member [%d].x \n", i);
            flag = 0;
        }
        if (!compare_float(seq_particles[i].y, par_particles[i].y)){
            printf("Difference in member [%d].y \n", i);
            flag = 0;
        }
        if (!compare_float(seq_particles[i].z, par_particles[i].z)){
            printf("Difference in member [%d].z \n", i);
            flag = 0;
        }
    }

    if (flag) printf("CORRETTO-L'output della computazione parallela è uguale all'output della compuzione sequenziale\n");
    else printf("NON CORRETTO-L'output della computazione parallela è diverso all'output della compuzione sequenziale\n");

    /* Infine chiudo i files per rilasciare le risorse */
    fclose(seq_file);
    fclose(par_file);
    free(seq_particles);
    free(par_particles);

    return 0;
}