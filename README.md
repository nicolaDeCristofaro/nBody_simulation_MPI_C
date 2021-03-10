
| N-Body Simulation | Nicola De Cristofaro |
| --- | --- |

# Problem Description
In physics, the n-body problem consists in predicting the individual movements of a group of celestial objects that interact with each other in a gravitational way. The resolution of this problem was motivated by the desire to understand the movements of the Sun, Moon, planets and visible stars.

<div style="text-align:center">
    <img alt="nBodySimulation1" src="./images/nbody1.gif" height="232" width="328"/>
    <img alt="nBodySimulation2" src="./images/nbody2.gif" height="232" width="328"/>
</div>

# Solution Description
To find the solution to this problem, it is possible to simulate the behavior of particles, each having an initial mass, position and velocity. The simulation will allow you to calculate the position and speed of each particle at the end of a given time interval.

The proposed solution has a quadratic complexity with respect to the size of the input (the number of particles). There is a type of simulation, called Barnes - Hut, which is more efficient as it is able to execute with order O (n log n) through approximations but it is more complex to implement therefore it was preferred to consider the simpler quadratic complexity solution and concentrate on concepts of parallel programming.

# Software used

- <b> Message Passing Interface (MPI) </b>:  A library used to allow several different processors on a cluster to communicate with each other. In other words, the goal of the Message Passing Interface is to provide a widely used standard for writing message-passing programs.

- <b> C language</b>

- <b> Amazon Web Services (AWS)</b>: to test the parallel program on a cluster of virtual machine create on AWS Cloud environment. 

# Solution Structure

## Initialization
For the initialization phase, the technique of creating the particles in a file and initializing them using a deterministic algorithm to randomize the values of the particles was adopted. In this way all the processors at the first iteration of the computation start reading from this file *(particles.txt)*.

To create and initialize the particles, execute the following commands:

1. ```bash
    gcc -o initialization particles_production.c
    ```
2. ```bash
    ./initialization [number of particles] (example -> ./initialization 1000)
    ```

## MPI Initialization, variables and memory allocation
The first part of the program deals with the definition of variables useful for computation, variables to take execution times and finally variables for the use of MPI. In addition to the operations always present such as MPI_Init, MPI_Comm_size to know the number of processors that are running or MPI_Comm_rank to know the current processor rank, another operation was made to create a derived data type in order to allow communication between processors of the <b>Particle</b> data type.

```c
    MPI_Type_contiguous(7, MPI_FLOAT, &particle_type);
    MPI_Type_commit(&particle_type);
```

where the struct Particle is the following:

```c
    typedef struct {
        float mass;
        float x, y, z;
        float vx, vy, vz;
    } Particle;
```

## Equal distribution of workload between processors
The function is then performed to distribute the work equally among the processors. After executing this function in the * dim_portions * array each index is associated with the rank of a processor and its contents represent the size of the portion of the array on which that processor will have to compute. Furthermore, the array of displacements is set in which each index represents the start_offset of a processor.

```c
    void compute_equal_workload_for_each_task(int *dim_portions, int *displs, int arraysize, int numtasks){

        for(int i=0; i<numtasks;i++){
            dim_portions[i] = (arraysize / numtasks) +
                            ((i < (arraysize % numtasks)) ? 1 : 0);
        }

        int offset = 0;
        for(int i=0;i<numtasks;i++){
            displs[i] = offset;
            offset += dim_portions[i];
        }
    }
```

## Computing and communication between processors

- The simulation takes place for a certain number of iterations I set to 10. Furthermore we consider the processor with rank 0 the MASTER processor while the others are considered SLAVES.

- The operations carried out in each iteration are the following:

1. In the first iteration only, all processors read the initial state of the particles from the file. If, on the other hand, it is not the first iteration, it means that the MASTER processor has the particle array updated after the computation of the previous iteration, so it broadcasts it to all the other processors for the new computation.

```c
        if(iteration == 0){
            //First iteration: all the procecssors can read the initial state of the particles form file
            FILE *fileRead = fopen("particles.txt", "r");
            if (fileRead == NULL){
                /* Impossible to open the file */
                printf("\nImpossible to open the file.\n");
                exit(EXIT_FAILURE);
            }

            fread(particles, sizeof(Particle) * num_particles, 1, fileRead);
            fclose(fileRead);
        }else{
            //MASTER processor has the output particle array of the previous computation so it broadcasts
            MPI_Bcast(particles, num_particles, particle_type, MASTER, MPI_COMM_WORLD);
        }
```

2. Each processor performs the computation on its portion of the array whose size was calculated fairly in the previous step by calling the * bodyForce * function which allows to calculate the new position and velocity values of each particle. ** It should be emphasized that the MASTER processor does not only manage the communication between processors but also performs the computation on its portion. **

```c
    bodyForce(particles, displ[myrank], dt, dim_portions[myrank], num_particles );
```

- as we can see from the function code, each processor performs the computation on a subset of particles whose values are computed based on the force exerted by all the other particles.

```c
 /*Function that performs computation on a specific portion of the workload */
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

    //Integration of the positions of my portion
    for(int i = 0; i < dim_portion; i++) {
        all_particles[startOffsetPortion + i].x += all_particles[startOffsetPortion + i].vx * dt;
        all_particles[startOffsetPortion + i].y += all_particles[startOffsetPortion + i].vy * dt;
        all_particles[startOffsetPortion + i].z += all_particles[startOffsetPortion + i].vz * dt;
    }
}
```

3. Each processor sends its computed portion to the MASTER so, after this gathering operation, the MASTER processor has the particle array complete and computed for a certain iteration.

```c
    MPI_Gatherv(particles + displ[myrank], dim_portions[myrank], particle_type, gathered_particles, dim_portions, displ, particle_type,MASTER, MPI_COMM_WORLD);
```

4. The input of the next iteration must be the particle array computed in the current iteration so since the MASTER processor owns the particle array computed in * gathered_particles * a swap is performed.

```c
    if(myrank == MASTER) particles = gathered_particles;
```

5. Finally, for each iteration the execution time is taken and the MASTER processor writes on stdout the progress of the computation and the time taken for that iteration.

```c
    MPI_Barrier(MPI_COMM_WORLD);
    iterStart = MPI_Wtime();
    ...
    ...
    MPI_Barrier(MPI_COMM_WORLD);  
    iterEnd = MPI_Wtime();
    if(myrank == MASTER) printf("Iteration %d of %d completed in %f seconds\n", iteration+1, I, (iterEnd-iterStart));
```

## Finalization and deallocation
In this last phase the computation is completed for all the iterations, then all the memory previously allocated is deallocated, MPI is finalized with MPI_Finalize () and the MASTER processor writes the average execution time of an iteration to stdout and the total execution time, also writing the final state of the particles on a file for a possible correctness test done later.

## Compilazione ed esecuzione

#### Compilazione 
```bash
    mpicc -o parallel parallel_nBody.c -lm
```

#### Esecuzione 
```bash
    mpirun -np [numero di processori] parallel [numero di particelle]
    esempio -> mpirun -np 4 parallel 1000
```
**ASSICURARSI CHE IL FILE DA CUI VENGONO LETTE LE PARTICELLE SIA STATO CREATO**

# Correttezza

- Per la correttezza di un programma in parallelo è necessario che l'esecuzione con P processori o con un solo processore, con lo stesso input produca lo stesso output. Per verificare questa condizione è stato fatto in modo che l'output dell'esecuzione parallela con più processori sia stato scritto su file *(parallel_output.txt)*.

- Eseguiamo ora la versione sequenziale del programma che può essere eseguita sia utilizzando la versione parallela su 1 processore, sia utilizzando la versione sequenziale che evita anche l'inizializzazione e finalizzazione di MPI.

1. ```bash
    gcc -o sequential sequential_nBody.c -lm
    ```

2. ```bash
    ./sequential [numero di particelle] 
    esempio -> ./sequential 1000
    ```
    **OPPURE**

 1. ```bash
    mpicc -o parallel parallel_nBody.c -lm
    ```

2. ```bash
    mpirun -np 1 parallel [numero di particelle]
    esempio -> mpirun -np 1 parallel 1000
    ```

- Per eseguire il test di correttezza utilizziamo la versione sequenziale del programma. L'output della versione sequenziale verrà scritto sul file *(sequential_output.txt)*. E' stato poi realizzato un programma *(output_correctness.c)* che mette a confronto il contenuto del file con l'output sequenziale e del file con l'output parallelo per verificarne la correttezza.

- Dopo aver eseguito il programma sia nella versione parallela che quella sequenziale come indicato precedentemente è possibile eseguire il test di correttezza nel seguente modo:

1. ```bash
    gcc -o correctness output_correctness.c
    ```

2. ```bash
    ./correctness [numero di particelle] (esempio -> ./correctness 1000)
    ```

- Da sottolineare che per quanto riguarda i confronti tra gli attributi delle particelle per stabilirne la correttezza è stata utilizzata una funzione per confrontare i valori float. Questo perchè la matematica in virgola mobile non è esatta. Valori semplici come 0,2 non possono essere rappresentati con precisione usando numeri binari in virgola mobile e questa precisione limitata dei numeri in virgola mobile significa che lievi variazioni nell'ordine delle operazioni possono cambiare il risultato. Quindi, dato che questi valori vengono memorizzati con precisioni diverse, i risultati possono differire. Di conseguenza se si eseguono vari calcoli e quindi si confrontano i risultati con un valore atteso, è altamente improbabile che si ottenga esattamente il risultato desiderato. Per questo motivo con la funzione realizzata possiamo esprimere il concetto che due valori abbastanza vicini tra loro possono essere considerati uguali.

- In pratica, due valori di tipo float sono considerati uguali se la loro differenza rientra in un certo limite o valore epsilon.

```c
//confronta due float e ritorna 1 se sono uguali 0 altrimenti
int compare_float(float f1, float f2){

    float precision = 0.01f;
    if (((f1 - precision) < f2) &&
        ((f1 + precision) > f2))
        return 1;
    else
        return 0;
}
```

- **Visualizzazione grafica della correttezza della soluzione parallela su input piccolo (5)**

![](./images/correctness_5.png)

# Problem evaluation and Benchmarks

Per la valutazione delle prestazioni della soluzione proposta sono state create delle versioni dei programmi **(sequenziale: sequential_nBody_benchmarking.c - parallela: parallel_nBody_benchmarking.c)** leggermente revisionate poichè per effettuare un migliore benchmarking è stata eliminata la scrittura dell'output finale su file per focalizzarci sul tempo computazione.

Possiamo ora procedere con la descrizione dei risultati dati dalla misurazione della scalabilità dell'applicazione. Ci sono due modi di base per misurare la performance parallela di un'applicazione: **strong e weak scaling**.

## Strong Scaling

In questo tipo di misurazione la taglia del problema (il numero di particelle) resta fissata ma il numero di processori aumenta. 

![](./benchmarking_screenshots/instances_32.png)

**Per il testi di strong scaling il range di aumento del numero dei processori è stato fissato da 1 a 32 poichè 32 è il massimo numero di core che è possibile usare con un account AWS Educate**

![](./benchmarking_screenshots/cpu_32_limit.png)

|   P	|   N	|   Avg Iteration Time (seconds)	|   Total Computation Time (seconds)	|   Strong Scaling Efficiency (%)  |
|:-:	|:-:	|:-:	|:-:	|:-:	|
|   1	|  100 000	|   220.005443	|   2200.054430 ≅ 37 min	| 100   |
|   2	|  100 000 	|   109.134224 	|   1091.342241 	|   100 |
|   3	|  100 000	|   73.359361	|   733.593608 	|   99  |
|   4	|  100 000 	|   54.819589	|   548.195891  |   100 |
|   5	|  100 000 	|   43.785572	|   437.855721	|   100 |
|   6	|  100 000 	|   36.821113	|   368.211132 	|   99  |
|   7	|  100 000 	|   31.567646	|   315.676457 	|   99  |
|   8	|  100 000 	|   27.412049	|   274.120492  |   100 |
|   9	|  100 000 	|   24.375618	|   243.756177 	|   100 |
|   10	|  100 000 	|   21.972933	|   219.729333 	|   100 |
|   11	|  100 000 	|   20.238017	|   202.380169 	|   98  |
|   12	|  100 000 	|   18.540094	|   185.400640 	|   98  |
|   13	|  100 000 	|   17.049654	|   170.496536 	|   99  |
|   14	|  100 000 	|   15.787466	|   157.874664 	|   99  |
|   15	|  100 000 	|   14.750758	|   147.507577 	|   99  |
|   16	|  100 000 	|   13.837183	|   138.371826 	|   99  |
|   17	|  100 000 	|   12.937212	|   129.372124	|   100 |
|   18	|  100 000 	|   12.210964	|   122.109644 	|   100 |
|   19	|  100 000 	|   11.575181	|   115.751814  |   100 |
|   20	|  100 000 	|   10.975263	|   109.752627 	|   100 |
|   21	|  100 000 	|   10.547308	|   105.473079  |   99  |
|   22	|  100 000 	|   9.991014	|   99.910142   |   100 |
|   23	|  100 000 	|   9.599492	|   95.994917   |   99  |
|   24	|  100 000 	|   9.162451	|   91.624507   |   100 |
|   25	|  100 000 	|   8.797949	|   87.979485   |   100 |
|   26	|  100 000 	|   8.514978	|   85.149779   |   99  |
|   27	|  100 000 	|   8.152149	|   81.521492   |   99  |
|   28	|  100 000 	|   7.863964	|   78.639643   |   99  |
|   29	|  100 000 	|   7.639826	|   76.398257   |   99  |
|   30	|  100 000 	|   7.336128	|   73.361283   |   99  |
|   31	|  100 000 	|   7.095153	|   70.951534   |   100 |
|   32	|  100 000 	|   6.880448	|   68.804478   |   99  |

![](./benchmarking_screenshots/strong_scaling_screenshots/strong_scaling_graph.png)

Dalle misurazioni effettuate per verificare la **strong scalability** abbiamo notato come su una stessa taglia di input, l'aggiunta di processori ha migliorato il tempo di esecuzione, ma ovviamente per ogni processore aggiunto il miglioramento non è stato costante poichè più processori partecipavano alla computazione, maggiore era l'overhead prodotto dalla comunicazione di questi. Nel nostro test fino a 32 processori il tempo di esecuzione è sempre diminuito ma più il numero di processori aumentava più il miglioramento delle prestazioni diminuiva. Se fossimo andati avanti, aggiungendo altri processori saremmo arrivati ad un punto in cui il tempo di esecuzione non diminuiva più, ma cominciava ad aumentare poichè l'overhead prodotto era maggiore del miglioramento di prestazioni.

## Weak Scaling

In questo caso la taglia del problema aumenta con l'aumentare del numero di processori, facendo in modo che il workload sia sempre equamente distribuito tra i processori.

|   P	|   N	|   Avg Iteration Time (seconds)	|   Total Computation Time (seconds)	|   Weak Scaling Efficiency (%)  |
|:-:	|:-:	|:-:	|:-:	|:-:	|
|   1	|  10 000	|   1.961227	|   19.612270	|   100 |
|   2	|  20 000 	|   3.921892	|   39.218918	|   50  |
|   3	|  30 000	|   5.893506	|   58.935057	|   33  |
|   4	|  40 000 	|   7.872926	|   78.729259	|   24  |
|   5	|  50 000 	|   10.999243	|   109.992425	|   17  |
|   6	|  60 000 	|   13.164654	|   131.645642	|   14  |
|   7	|  70 000 	|   15.374346	|   153.743464	|   12  |
|   8	|  80 000 	|   17.566746	|   175.667462	|   11  |
|   9	|  90 000 	|   19.766140	|   197.661402	|   9   |
|   10	|  100 000 	|   21.990926	|   219.909259	|   8   |
|   11	|  110 000 	|   24.507350	|   245.073497	|   8   |
|   12	|  120 000 	|   26.600296	|   266.022962	|   7   |
|   13	|  130 000 	|   28.731199	|   287.311988	|   6   |
|   14	|  140 000 	|   30.937833	|   309.378334	|   6   |
|   15	|  150 000 	|   33.151154	|   331.511537	|   5   |
|   16	|  160 000 	|   35.357112	|   353.571116	|   5   |

![](./benchmarking_screenshots/weak_scaling_screenshots/weak_scaling_graphic.png)

Il grafico ideale della performance di **weak scalability** sarebbe una linea retta poichè la taglia dell'input è aumentata in proporzione all'aumento del numero di processori, quindi essendo il workload sempre equamente distrubuito il tempo di computazione dovrebbe essere sempre lo stesso. Purtroppo questo non accade poichè, come abbiamo già detto, per ogni processore aggiunto viene prodotta una quantità maggiore di overhead, dovuto principalmente alla comunicazione tra i processori. Dal nostro esperimento abbiamo avuto comunque buoni risultati poichè con l'aumentare dei processori e della taglia dell'input (10.000 particelle in più per ogni processore aggiunto alla computazione) il tempo di esecuzione è aumentato in modo minimo e costante ad ogni step, rimanendo relativamente vicino all'ideale.

## Speedup ed efficienza
Lo speedup è un'altra metrica di misurazione delle performance che rappresenta il miglioramento delle prestazioni di un programma dovuto a un'esecuzione parallela rispetto a una sequenziale.

Di solito il meglio in cui possiamo sperare è di dividere il lavoro tra i processori senza aggiungere ulteriore tempo di lavoro extra. In una situazione del genere se lanciamo il nostro programma con P processori allora il programma eseguirà P volte più velocemente del programma sequenziale. Se questo accade diciamo che il programma parallelo ha **speedup lineare**.

In pratica, difficilmente riusciamo ad ottenere questo tipo di speedup poichè l'uso di più processori introduce nella maggior parte dei casi un tempo di "overhead" ovvero un tempo extra impiegato principalmente per operazioni di comunicazione tra i processori.

L'efficienza invece è definita come il rapporto tra speedup e numero di processori e misura la frazione di tempo durante il quale un processore viene utilmente utilizzato.

Di seguito i valori relativi allo speedup e all'efficienza calcolati nel seguente modo:

**Speedup**

$$ S(P,N) = \frac{T(1,N)}{T(P,N)} $$

**Efficiency**

$$ E = \frac{S}{P} = \frac{\frac{T(1,N)}{T(P,N)}}{P} = \frac{T(1,N)}{p * T(P,N)} $$

dove :
+ **P = numero di processori usati**
+ **N = taglia dell'input**


|   P	|   N	|   Speedup	| Efficiency |
|:-:	|:-:	|:-:	|:-:	|
|   1	|  100 000	|  1.0 | 1.00 |
|   2	|  100 000 	|  2.0 | 1.00 |
|   3	|  100 000	|  3.0 | 1.00 |
|   4	|  100 000 	|  4.0 | 1.00 |
|   5	|  100 000 	|  5.0 | 1.00 |
|   6	|  100 000 	|  5.9 | 1.00 |
|   7	|  100 000 	|  6.9 | 1.00 |
|   8	|  100 000 	|  8.0 | 1.00 |
|   9	|  100 000 	|  9.0 | 1.00 |
|   10	|  100 000 	|  10.0 | 1.00|
|   11	|  100 000 	|  10.9 | 0.99|
|   12	|  100 000 	|  11.9 | 0.99|
|   13	|  100 000 	|  12.9 | 0.99|
|   14	|  100 000 	|  13.9 | 0.99|
|   15	|  100 000 	|  14.9 | 1.00|
|   16	|  100 000 	|  15.9 | 0.99|
|   17	|  100 000 	|  17.0 | 0.99|
|   18	|  100 000 	|  18.0 | 1.00|
|   19	|  100 000 	|  19.0 | 1.00|
|   20	|  100 000 	|  20.0 | 1.00|
|   21	|  100 000 	|  20.9 | 1.00|
|   22	|  100 000 	|  22.0 | 0.99|
|   23	|  100 000 	|  22.9 | 1.00|
|   24	|  100 000 	|  24.0 | 1.00|
|   25	|  100 000 	|  25.0 | 1.00|
|   26	|  100 000 	|  25.8 | 0.99|
|   27	|  100 000 	|  27.0 | 1.00|
|   28	|  100 000 	|  28.0 | 1.00|
|   29	|  100 000 	|  28.8 | 0.9|
|   30	|  100 000 	|  30.0 | 1.00|
|   31	|  100 000 	|  31.0 | 1.00|
|   32	|  100 000 	|  32.0 | 1.00|

![](./benchmarking_screenshots/speedup_32.png)

![](./benchmarking_screenshots/efficiency.png)

Abbiamo notato come, sorprendentemente, sia lo **speedup** che l'**efficienza** sono sempre stati molto vicini all'ideale, questo perchè sono calcolati in funzione di quante operazioni riusciamo ad eseguire in parallelo e nel nostro caso abbiamo avuto ottimi risultati poichè, tranne per alcune istruzioni in più eseguite dal processore MASTER, tutte le istruzione del programma sono state eseguite in parallelo.

Inoltre la comunicazione tra i processori è effettuata interamente utilizzando funzioni di comunicazione collettiva della libreria MPI e questo è un altro fattore che migliora le prestazioni poichè fa in modo che qualsiasi nodo con le informazioni ricevute partecipi all'invio delle informazioni ad altri nodi.


# Conclusioni

Dalle considerazione effettuate sulla strong scalability, weak scalability e sullo speedup ed efficienza possiamo concludere che abbiamo raggiunto buoni risultati e notato evidenti miglioramenti con l'aggiunta di altri processori. Infatti essendo un problema di complessità quadratica il tempo di esecuzione sequenziale era davvero elevato (su taglia di input 100 000 il tempo totale di esecuzione era circa 37 minuti). Impiegando 32 processori in parallelo il tempo totale di computazione è diminuito fino a circa 1 minuto, notando così la potenza di una computazione parallela.

<script type="text/javascript" src="http://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML"></script>
<script type="text/x-mathjax-config">
    MathJax.Hub.Config({ tex2jax: {inlineMath: [['$', '$']]}, messageStyle: "none" });
</script>



