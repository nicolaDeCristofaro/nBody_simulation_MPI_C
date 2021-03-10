
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

## Compilation and execution

#### Compilation 
```bash
    mpicc -o parallel parallel_nBody.c -lm
```

#### Execution 
```bash
    mpirun -np [number of processors] parallel [number of particles]
    example -> mpirun -np 4 parallel 1000
```
** MAKE SURE THE FILE FROM WHICH PARTICLES ARE READ HAS BEEN CREATED **

# Correctness

- For the correctness of a parallel program it is necessary that the execution with P processors or with a single processor, with the same input produces the same output. To verify this condition, the output of the parallel multi-processor execution was written to file * (parallel_output.txt) *.

- We now run the sequential version of the program which can be run either using the parallel version on 1 processor, or using the sequential version which also avoids the initialization and finalization of MPI.

1. ```bash
    gcc -o sequential sequential_nBody.c -lm
    ```

2. ```bash
    ./sequential [number of particles] 
    example -> ./sequential 1000
    ```
    **OR**

 1. ```bash
    mpicc -o parallel parallel_nBody.c -lm
    ```

2. ```bash
    mpirun -np 1 parallel [number of particles] 
    example -> mpirun -np 1 parallel 1000
    ```

- To perform the correctness test we use the sequential version of the program. The output of the sequential version will be written to the file * (sequential_output.txt) *. A program * (output_correctness.c) * was then created that compares the contents of the file with the sequential output and of the file with the parallel output to verify its correctness.

- After running the program in both the parallel and sequential versions as indicated above, it is possible to perform the correctness test as follows:

1. ```bash
    gcc -o correctness output_correctness.c
    ```

2. ```bash
    ./correctness [number of particles] (example -> ./correctness 1000)
    ```

- It should be emphasized that as regards the comparisons between the attributes of the particles to establish their correctness, a special function was used to compare the float values. This is because floating point math is not exact. Simple values such as 0.2 cannot be accurately represented using floating point binary numbers, and this limited precision of floating point numbers means that slight variations in the order of operations can change the result. Therefore, since these values are stored with different accuracies, the results may differ. Consequently, if you perform various calculations and then compare the results with an expected value, it is highly unlikely that you will get exactly the desired result. For this reason, with the realized function we can express the concept that two values close enough to each other can be considered equal.

- In practice, two float values are considered equal if their difference falls within a certain epsilon limit or value.
- 
```c
// compare two floats and return 1 if they are equal to 0 otherwise
int compare_float(float f1, float f2){

    float precision = 0.01f;
    if (((f1 - precision) < f2) &&
        ((f1 + precision) > f2))
        return 1;
    else
        return 0;
}
```

- ** Graphic display of the correctness of the parallel solution on small input (5) **
- 
![](./images/correctness_5.png)

# Problem evaluation and Benchmarks

To evaluate the performance of the proposed solution, some versions of the programs have been created ** (sequential: sequential_nBody_benchmarking.c - parallel: parallel_nBody_benchmarking.c) ** slightly revised, since to perform a better benchmarking the writing of the final output has been eliminated on file to focus on computation time.

We can now proceed with the description of the results given by measuring the scalability of the application. There are two basic ways to measure the parallel performance of an application: ** strong and weak scaling **.

## Strong Scaling

In this type of measurement the size of the problem (the number of particles) remains fixed but the number of processors increases.

![](./benchmarking_screenshots/instances_32.png)

** For the strong scaling test the range of increase in the number of processors has been fixed from 1 to 32 as 32 is the maximum number of cores that can be used with an AWS Educate account **

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

From the measurements carried out to verify the ** strong scalability ** we noticed how on the same input size, the addition of processors improved the execution time, but obviously for each processor added the improvement was not constant as more processors participated to computation, the greater was the overhead produced by the communication of these. In our test up to 32 processors the execution time always decreased but the more the number of processors increased the more the performance improvement decreased. If we had moved on, adding more processors we would have reached a point where the running time no longer decreased, but began to increase as the overhead produced was greater than the performance improvement.

## Weak Scaling

In this case the size of the problem increases as the number of processors increases, ensuring that the workload is always equally distributed among the processors.

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

The ideal graph of the ** weak scalability ** performance would be a straight line since the input size has increased in proportion to the increase in the number of processors, so since the workload is always equally distributed, the computation time should always be the same . Unfortunately this does not happen because, as we have already said, a greater amount of overhead is produced for each added processor, mainly due to the communication between the processors. However, from our experiment we had good results because with the increase of the processors and the input size (10,000 more particles for each processor added to the computation) the execution time has increased in a minimal and constant way at each step, remaining relatively close to ideal.

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



