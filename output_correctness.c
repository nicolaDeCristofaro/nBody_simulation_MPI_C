#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <errno.h>

typedef struct{
    float mass;
    float x, y, z;
    float vx, vy, vz;
} Particle;

/* Definzione di funzioni */
int compare_float(float f1, float f2);
int convertStringToInt(char *str);


/*Programma per valutare la correttezza dell'implementazione parallela:
Confronta l'output della versione sequenziale con l'output della soluzione parallela*/
int main(int argc, char *argv[]){

    FILE *seq_file = fopen("sequential_output.txt", "r");
    FILE *par_file = fopen("parallel_output.txt", "r");

    int num_particles = 1000;  //Numero delle particelle di DEFAULT se nessun parametro è fornito sulla command-line
    if(argc > 1){
        // E' stato fornito il parametro da linea di comando che indica il numero di particelle
        //Converti stringa a intero
        num_particles = convertStringToInt(argv[1]);
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

    /* Leggo lo stato finale delle particelle dalla computazione sequenziale e da quella parallela */
    int seq_particles_read = fread(seq_particles, sizeof(Particle) * num_particles, 1, seq_file);
    int par_particles_read = fread(par_particles, sizeof(Particle) * num_particles, 1, par_file);

    /* Controllo errori nel numero di particelle lette */
    if( seq_particles_read == 0){
        /*il numero di particelle da leggere è maggiore del numero di particelle nel file*/
        printf("ERROR: Il numero di particelle da leggere è maggiore del numero di particelle nel file (sequential_output.txt)\n");
        exit(EXIT_FAILURE);
    }
    if( par_particles_read == 0){
        /*il numero di particelle da leggere è maggiore del numero di particelle nel file*/
        printf("ERROR: Il numero di particelle da leggere è maggiore del numero di particelle nel file (arallel_output.txt)\n");
        exit(EXIT_FAILURE);
    }

    /*Verifico se lo stato finale delle particelle dopo la computazione sequenziale 
    è uguale a quello delle particelle dopo la computazione parallela*/
    int flag = 1;   //flag rimarrà uguale a 1 se i due autput sono uguali altrimenti sarà 0
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
        if (!compare_float(seq_particles[i].vx, par_particles[i].vx)){
            printf("Difference in member [%d].vx \n", i);
            flag = 0;
        }
        if (!compare_float(seq_particles[i].vy, par_particles[i].vy)){
            printf("Difference in member [%d].vy \n", i);
            flag = 0;
        }
        if (!compare_float(seq_particles[i].vz, par_particles[i].vz)){
            printf("Difference in member [%d].vz \n", i);
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

/* Confronta due float e ritorna 1 se sono uguali 0 altrimenti*/
int compare_float(float f1, float f2){

    float precision = 0.01f;
    if (((f1 - precision) < f2) &&
        ((f1 + precision) > f2))
        return 1;
    else
        return 0;
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