#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <string.h>

#define MASTER 0            // Rank of the Master processor  
#define I 10                // Number of iterations
#define SOFTENING 1e-9f     // Infinitely large value     

typedef struct {
    float mass;
    float x, y, z;
	float vx, vy, vz;
} Particle;


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

/*Funzione che esegue computazione su una specifica porzione di workload */
void bodyForce(Particle *particles_portion, Particle *particles_portion_acting, float dt, int *dim_portions, int round, int myrank, int numtasks) {
    for (int i = 0; i < dim_portions[myrank]; i++) { 
        float Fx = 0.0f; float Fy = 0.0f; float Fz = 0.0f;

        for (int j = 0; j < dim_portions[(myrank + round) % numtasks]; j++) {
            float dx = particles_portion_acting[j].x - particles_portion[i].x;
            float dy = particles_portion_acting[j].y - particles_portion[i].y;
            float dz = particles_portion_acting[j].z - particles_portion[i].z;
            float distSqr = dx*dx + dy*dy + dz*dz + SOFTENING;
            float invDist = 1.0f / sqrtf(distSqr);
            float invDist3 = invDist * invDist * invDist;

            Fx += dx * invDist3; Fy += dy * invDist3; Fz += dz * invDist3;
        }

        particles_portion[i].vx += dt*Fx; 
        particles_portion[i].vy += dt*Fy; 
        particles_portion[i].vz += dt*Fz;
    }
}

int main(int argc, char* argv[]){

    MPI_Datatype particle_type;             // MPI datatype per comunicare il tipo di dato "Particle"  
    int numtasks;                           // Numero di processori usati
    int myrank;                             // Rank del processo corrente
    double start, end, iterStart, iterEnd;  // Variabili usate per la misurazione del tempo di esecuzione totale e di ogni iterazione

    int *dim_portions;                      // Dimensione della porzione di workload per ogni processo
    int *displ;                             // Offset di partenza della porzione  di workload per ogni processo
    Particle *my_portion;                   // Porzione di particelle di un processo
    int left_rank,right_rank;               // Rank dei vicini


    int num_particles = 1000;  //Numero delle particelle di DEFAULT se nessun parametro è fornito sulla command-line
    if(argc > 1){
        // E' stato fornito il parametro da linea di comando che indica il numero di particelle
        num_particles = atoi(argv[1]);
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

    /*** Il processo MASTER inizializza le particelle ***/
    Particle *particles = NULL;
    if(myrank == MASTER) {
	    particles = (Particle*) malloc(num_particles * sizeof(Particle));

        srand(myrank);
        randomizeParticles(particles,num_particles);

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


    /*** Distribuzione delle porzioni ai vari processi ***/
    my_portion = (Particle*) malloc(sizeof(Particle) * dim_portions[myrank]);
	MPI_Scatterv(particles, dim_portions, displ, particle_type,my_portion, dim_portions[myrank], particle_type,MASTER, MPI_COMM_WORLD);

    /*** Inizio della simulazione ***/
	Particle *portion_send;
	Particle *portion_recv;
    int portion_size = dim_portions[MASTER]; // Una porzione da inviare o da ricevere non può essere di dimensioni maggiori di quella del MASTER
    
    portion_send = (Particle*) malloc(sizeof(Particle) * portion_size);
	portion_recv = (Particle*) malloc(sizeof(Particle) * portion_size);

    const float dt = 0.01F; // time step
	for(int iteration = 0; iteration < I; iteration++) {

        MPI_Barrier(MPI_COMM_WORLD);  //Sincronizzo i processi prima di cominciare a prendere il tempo di esecuzione dell'iterazione
        iterStart = MPI_Wtime();	    
	  
      	/*** Computazione eseguita per calcolare la forza esercitata da tutte le particelle sulle particelle della mia porzione*/
	    for(int round = 0; round < numtasks; round++) {
		
            MPI_Request request[2];
            MPI_Status status[2];

            // Se siamo nel primo round la porzione du cui lavorare è quella locale
            if(round == 0) memcpy(portion_send, my_portion, dim_portions[myrank] * sizeof(Particle));

            //Cerco il rank dei miei vicini
            if ((left_rank=(myrank-1)) < 0) left_rank = numtasks-1;
            if ((right_rank=(myrank+1)) == numtasks) right_rank = 0; 

            // Invio asincrono ai miei vicini nel round corrente
            if(round < (numtasks - 1)) {
                //Invio al processo di sinistra 
                MPI_Isend(portion_send, portion_size, particle_type, left_rank,99, MPI_COMM_WORLD, &request[0]);
                //Ricevo dal processo di destra
                MPI_Irecv(portion_recv, portion_size, particle_type, right_rank,99, MPI_COMM_WORLD, &request[1]);
            }

            //Computazione forza esercitata sulle particelle 		
            bodyForce(my_portion, portion_send, dt, dim_portions, round, myrank, numtasks );

            //Attendo il completamento delle send/receive asincrone
            if(round < (numtasks - 1)) {
                MPI_Waitall(2, request, status);
                //Copio la porzione ricevuta da destra per computare la forza esercitata al prossimo round
                memcpy(portion_send, portion_recv,portion_size * sizeof(Particle));
            }
        }

        /*** Tutte le particelle conoscono la forza esercitata dalle altre particelle quindi 
          è possibile calcolare le nuove posizioni delle particelle ***/
        for(int i = 0; i < dim_portions[myrank]; i++) { 
            my_portion[i].x += my_portion[i].vx*dt;
            my_portion[i].y += my_portion[i].vy*dt;
            my_portion[i].z += my_portion[i].vz*dt;
        }

        iterEnd = MPI_Wtime();
        if(myrank == MASTER) printf("Iterazione %d di %d completata in %f seconds\n", iteration, I, (iterEnd-iterStart));
    }

    /*** Gathering della porzione computata da ogni processo ***/
	Particle *gathered_particles = NULL;
	if(myrank == MASTER) gathered_particles = (Particle*) malloc(sizeof(Particle) * num_particles);
	
    MPI_Gatherv(my_portion, dim_portions[myrank], particle_type, gathered_particles, dim_portions, displ, particle_type,MASTER, MPI_COMM_WORLD);
	
    free(portion_send);
	free(portion_recv);
	free(my_portion);
    free(dim_portions);
    free(displ);

    MPI_Barrier(MPI_COMM_WORLD);     // tutti i processi hanno terminato 
    end = MPI_Wtime();               // Prendo il tempo di fine esecuzione
    MPI_Finalize();                  

    if(myrank == MASTER) {
        double totalTime = end-start;
        double avgTime = totalTime / (double)(I); 
        printf("\nAvg iteration time: %f seconds\n", avgTime);
        printf("Total time: %f\n", totalTime);

        /* TEST: decommenta per scrivere su stdout lo stato finale delle particelle dopo la computazione
        printf("\nOUTPUT\n");
        for (int i = 0  ; i < num_particles; i++) { 
            printf("[%d].x = %f\t", i, gathered_particles[i].x);
            printf("[%d].y = %f\t", i, gathered_particles[i].y);
            printf("[%d].z = %f\t", i, gathered_particles[i].z);
            printf("[%d].vx = %f\t", i, gathered_particles[i].vx);
            printf("[%d].vy = %f\t", i, gathered_particles[i].vy);
            printf("[%d].vz = %f\t", i, gathered_particles[i].vz);
            printf("\n");
        }*/

        /*Scrivo l'output su file per poi poterne valutare la correttenza confrontando con l'output sequenziale*/
        FILE * file= fopen("parallel_output.txt", "wb");
        if (file != NULL) {
            fwrite(gathered_particles , sizeof(Particle) * num_particles, 1, file);
            fclose(file);
        }
	    
        free(particles);
	    free(gathered_particles);	
	}
    
    return 0;
}
