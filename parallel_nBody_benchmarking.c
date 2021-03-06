#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <errno.h>
#include <mpi.h>

#define MASTER 0            // Rank del processore MASTER 
#define I 10                // Numero di iterazioni
#define SOFTENING 1e-9f     // Valore infinitamente grandeutilizzato nella computazione  

typedef struct {
    float mass;
    float x, y, z;
	float vx, vy, vz;
} Particle;


//Definizione delle funzioni
void compute_equal_workload_for_each_task(int *dim_portions, int *displs, int arraysize, int numtasks);
void bodyForce(Particle *all_particles, int startOffsetPortion, float dt, int dim_portion, int num_particles);
int convertStringToInt(char *str);


int main(int argc, char* argv[]){

    MPI_Datatype particle_type;             // MPI datatype per comunicare il tipo di dato "Particle"  
    int numtasks;                           // Numero di processori usati
    int myrank;                             // Rank del processo corrente
    double start, end, iterStart, iterEnd;  // Variabili usate per la misurazione del tempo di esecuzione totale e di ogni iterazione

    int *dim_portions;                      // Dimensione della porzione di workload per ogni processo
    int *displ;                             // Offset di partenza della porzione  di workload per ogni processo
    Particle *my_portion;                   // Porzione di particelle di un processo

    int num_particles = 1000;  //Numero delle particelle di DEFAULT se nessun parametro è fornito sulla command-line
    if(argc > 1){
        // E' stato fornito il parametro da linea di comando che indica il numero di particelle
        //Converti stringa a intero
        num_particles = convertStringToInt(argv[1]);
    }
        
    /*** Inizialliza MPI ***/
    MPI_Init(&argc, &argv);

    /*** Creazione del tipo di dato MPI per comunicare il tipo di dato "Particle" ***/
    MPI_Type_contiguous(7, MPI_FLOAT, &particle_type);
    MPI_Type_commit(&particle_type);
    
    /*** Ottengo il numero di processori usati e il rank del processo corrente ***/
    MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Barrier(MPI_COMM_WORLD);    /* tutti i processi sono inizializzati */
    start = MPI_Wtime();            /* prendo il tempo di inizio dell'esecuzione */

    /*** Calcolo come le particelle sono associate equamente ai processi ***/
    dim_portions = (int*) malloc(sizeof(int)*numtasks);
    displ = (int*) malloc(sizeof(int)*numtasks);
    compute_equal_workload_for_each_task(dim_portions, displ, num_particles, numtasks);

    const float dt = 0.01f; // time step
    
    Particle *particles = (Particle*) malloc(num_particles * sizeof(Particle));
    my_portion = (Particle*) malloc(sizeof(Particle) * dim_portions[myrank]);
    Particle *gathered_particles = NULL;
	if(myrank == MASTER) gathered_particles = (Particle*) malloc(sizeof(Particle) * num_particles);

    for(int iteration = 0; iteration < I; iteration++) {

        MPI_Barrier(MPI_COMM_WORLD);  //Sincronizzo i processi prima di cominciare a prendere il tempo di esecuzione dell'iterazione
        iterStart = MPI_Wtime();

        if(iteration == 0){
            //E' la prima iterazione quindi tutti i processori possono leggere lo stato iniziale delle particelle da file
            FILE *fileRead = fopen("particles.txt", "r");
            if (fileRead == NULL){
                /* Impossibile aprire il file */
                printf("\nImpossibile aprire il file.\n");
                exit(EXIT_FAILURE);
            }

            int particlesRead = fread(particles, sizeof(Particle) * num_particles, 1, fileRead);
            if( particlesRead == 0){
                /*il numero di particelle da leggere è maggiore del numero di particelle nel file*/
                printf("ERROR: Il numero di particelle da leggere è maggiore del numero di particelle nel file\n");
                exit(EXIT_FAILURE);
            }

            fclose(fileRead);
        }else{
            //il processore MASTER ha l'array di particelle output della computazine precedente quindi spedisce in broadcast
            MPI_Bcast(particles, num_particles, particle_type, MASTER, MPI_COMM_WORLD);
        }

        bodyForce(particles, displ[myrank], dt, dim_portions[myrank], num_particles );

        /*** Gathering della porzione computata da ogni processo ***/
        MPI_Gatherv(particles + displ[myrank], dim_portions[myrank], particle_type, gathered_particles, dim_portions, displ, particle_type,MASTER, MPI_COMM_WORLD);

        if(myrank == MASTER) particles = gathered_particles;

        MPI_Barrier(MPI_COMM_WORLD);  
        iterEnd = MPI_Wtime();
        if(myrank == MASTER) printf("Iterazione %d di %d completata in %f seconds\n", iteration+1, I, (iterEnd-iterStart));
    }

    MPI_Barrier(MPI_COMM_WORLD);     // tutti i processi hanno terminato 
    end = MPI_Wtime();               // Prendo il tempo di fine esecuzione
    MPI_Finalize();   

    if(myrank == MASTER) {
        double totalTime = end-start;
        double avgTime = totalTime / (double)(I); 
        printf("\nAvg iteration time: %f seconds\n", avgTime);
        printf("Total time: %f seconds\n", totalTime);
	}

    free(my_portion);
    free(dim_portions);
    free(displ);
    free(particles);
    
    return 0;

}

/*Distribuzione equa del lavoro tra i tasks*/
void compute_equal_workload_for_each_task(int *dim_portions, int *displs, int arraysize, int numtasks){

    for(int i=0; i<numtasks;i++){
        dim_portions[i] = (arraysize / numtasks) +
                        ((i < (arraysize % numtasks)) ? 1 : 0);
    }

    //imposto l'array dei displacements : ogni indice rappresenta lo start_offset di un task
    int offset = 0;
    for(int i=0;i<numtasks;i++){
        displs[i] = offset;
        offset += dim_portions[i];
    }

    /*dopo questa funzione nell'array dim_portions ogni indice è associato a un task 
    il valore associato a uno specifico indice rappresenta la dimensione della porzione di workload associata a quel task*/
}

/*Funzione che esegue computazione su una specifica porzione di workload */
void bodyForce(Particle *all_particles, int startOffsetPortion, float dt, int dim_portion, int num_particles) {
    for (int i = 0; i < dim_portion; i++) { 
        float Fx = 0.0f; float Fy = 0.0f; float Fz = 0.0f;

        for (int j = 0; j < num_particles; j++) {
            float dx = all_particles[j].x - all_particles[startOffsetPortion + i].x;
            float dy = all_particles[j].y - all_particles[startOffsetPortion + i].y;
            float dz = all_particles[j].z - all_particles[startOffsetPortion + i].z;
            float distSqr = dx*dx + dy*dy + dz*dz + SOFTENING;
            float invDist = 1.0f / sqrtf(distSqr);
            float invDist3 = invDist * invDist * invDist;

            Fx += dx * invDist3; Fy += dy * invDist3; Fz += dz * invDist3;
        }

        all_particles[startOffsetPortion + i].vx += dt * Fx; 
        all_particles[startOffsetPortion + i].vy += dt * Fy; 
        all_particles[startOffsetPortion + i].vz += dt * Fz;
    }

    //Integro le posizioni della mia porzione
    for(int i = 0; i < dim_portion; i++) {
        all_particles[startOffsetPortion + i].x += all_particles[startOffsetPortion + i].vx * dt;
        all_particles[startOffsetPortion + i].y += all_particles[startOffsetPortion + i].vy * dt;
        all_particles[startOffsetPortion + i].z += all_particles[startOffsetPortion + i].vz * dt;
    }
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