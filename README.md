| N-Body Simulation | Nicola De Cristofaro | Data di consegna |
| --- | --- | --- |

# Descrizione del problema
In fisica, il problema n-body consiste nel predire i  singoli movimenti di un gruppo di oggetti celesti che interagiscono tra loro in modo gravitazionale. La risoluzione di questo problema è stata motivata dal desiderio di comprendere i movimenti del Sole, della Luna, dei pianeti e delle stelle visibili.

<div style="text-align:center">
    <img alt="nBodySimulation1" src="./images/nbody1.gif" height="232" width="328"/>
    <img alt="nBodySimulation2" src="./images/nbody2.gif" height="232" width="328"/>
</div>

# Descrizione della soluzione
Per trovare la soluzione a questo problema è possibile simulare il comportamento delle particelle, ognuna avente una massa, una posizione e una velocità iniziale. La simulazione consentirà di calcolare la posizione e la velocità di ciascuna particella al termine di un intervallo di tempo determinato.

La soluzione proposta ha una complessità quadratica rispetto alla taglia dell'input (il numero di particelle). Esiste un tipo di simulazione, detta di Barnes – Hut che è più efficiente poichè tramite approssimazioni riesce ad eseguire con ordine O(n log n) ma è più complessa da implementare perciò si è preferito considerare la soluzione di complessità quadratica più semplice e concentrarsi sui concetti della programmazione parallela.

# Struttura della soluzione

### Inizializzazione
Per la fase di inizializzazione è stata adottata la tecnica di creare le particelle in un file e inizializzarle usando un algoritmo deterministico per randomizzare i valori delle particelle. In questo modo tutti i processori alla prima iterazione della computazione partono leggendo da questo file *(particles.txt)*.

Per la creazione e inizializzazione delle particelle eseguire i seguenti comandi:

1. ```bash
    gcc -o initialization particles_production.c
    ```
2. ```bash
    ./initialization [numero di particelle] (esempio -> ./initialization 1000)
    ```

### Inizializzazione di MPI, variabili e allocazione memoria
La prima parte del programma parte con la definizione di variabili utili alla computazione, variabili per prendere i tempi di esecuzione e infine variabili per l'uso di MPI. Oltre alle operazioni sempre presenti come MPI_Init, MPI_Comm_size per conoscere il numero di processori che stanno eseguendo o MPI_Comm_rank per conoscere il rank del processore corrente, è stata fatta un'altra operazione per creare un tipo di dato derivato al fine di permettere la comunicazione tra i processori del tipo di dato Particle.

```c
    MPI_Type_contiguous(7, MPI_FLOAT, &particle_type);
    MPI_Type_commit(&particle_type);
```

### Distribuziona equa del workload tra i processori
Viene quindi eseguita la funzione per distribuire equamente il lavoro tra i processori. Dopo l'esecuzione di questa funzione nell'array dim_portions ogni indice è associato al rank di un processore e il suo contenuto rappresenta la dimensione della porzione di array su cui quel processore dovrà computare. Inoltre viene impostato l'array dei displacements in cui ogni indice rappresenta lo start_offset di un processore.

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

### Computazione e comunicazione tra i processori

- La simulazione avviene per un certo numero di iterazioni I impostato a 10. Inoltre consideriamo il processore con rank 0 il processore MASTER mentre gli altri vengono considerati SLAVES.

- Le operazioni effettuate in ogni iterazione sono le seguenti:

1. Solo nella prima iterazione tutti i processori leggono lo stato iniziale delle particelle dal file. Se invece non è la prima iterazione significa che il processore MASTER ha l'array di particelle aggiornato dopo la computazione dell'iterazione precedente quindi lo spedisce in broadcast a tutti gli altri processori per la nuova computazione.

```c
    MPI_Bcast(particles, num_particles, particle_type, MASTER, MPI_COMM_WORLD);
```

2. Il processore MASTER spedisce inoltre a tutti i processorri la porzione di array la cui dimensione è stata calcolata in modo equo al passo precedente. **Da sottolineare che il processore MASTER non si occupa solo di gestire la comunicazione tra processori ma effettua anch'esso la computazione sulla sua porzione.**

```c
    MPI_Scatterv(particles, dim_portions, displ, particle_type,my_portion, dim_portions[myrank], particle_type,MASTER, MPI_COMM_WORLD);
```

3. Ogni processore effettua la computazione sulla sua porzione di array chiamando la funzione *bodyForce* che permette di calcolare i nuovi valori di posizione e velocità di ogni particella.

```c
    bodyForce(particles, my_portion, dt, dim_portions[myrank], num_particles );
```

4. Ogni processore spedisce al MASTER la proprio porzione computata quindi dopo questa operazione di gathering il processore MASTER ha l'array di particelle completo e computato per una certa iterazione.

```c
    MPI_Gatherv(my_portion, dim_portions[myrank], particle_type, gathered_particles, dim_portions, displ, particle_type,MASTER, MPI_COMM_WORLD);
```

5. L'input della successiva iterazione dovrà essere l'array di particelle computato nell'iterazione corrente quindi viene effettuato uno swap.

```c
    if(myrank == MASTER) particles = gathered_particles;
```

6. Infine per ogni iterazione viene preso il tempo di esecuzione e il processore MASTER provvede a scrivere su stdout lo stato di avanzamento della computazione e il tempo impiegato per quella iterazione.

```c
    MPI_Barrier(MPI_COMM_WORLD);
    iterStart = MPI_Wtime();
    ...
    ...
    MPI_Barrier(MPI_COMM_WORLD);  
    iterEnd = MPI_Wtime();
    if(myrank == MASTER) printf("Iterazione %d di %d completata in %f seconds\n", iteration, I, (iterEnd-iterStart));
```

### Finalizzazione e deallocazione
In quest'ultima fase la computazione è conclusa per tutte le I iterazioni quindi viene deallocata tutta la memoria precedentemente allocata, viene finalizzato MPI con MPI_Finalize() e il processore MASTER provvede a scrivere su stdout il tempo medio di esecuzione di un'iterazione e il tempo totale di esecuzione, scrivendo inoltre lo stato finale delle particelle su un file per un eventuale test di correttezza fatto in seguito.

### Compilazione ed esecuzione

#### Compilazione 
```bash
    mpicc -o parallel_nBody parallel_nBody_v4.0.c -lm
```

#### Esecuzione 
```bash
    mpirun -np [numero di processori] parallel_nBody [numero di particelle]
    esempio -> mpirun -np 4 parallel_nBody 1000
```

# Correttezza

- Per la correttezza di un programma in parallelo è necessario che l'esecuzione con P processori o con un solo processore, con lo stesso input produca lo stesso output. Per verificare questa condizione è stato fatto in modo che l'output dell'esecuzione parallela con più processori sia stato scritto su file *(parallel_output.txt)*.

- Eseguiamo ora la versione sequenziale del programma:

1. ```bash
    gcc -o sequential_nBody sequential_nBody.c -lm
    ```

2. ```bash
    ./sequential_nBody [numero di particelle] 
    esempio -> ./sequential_nBody 1000
    ```

- L'output della versione sequenziale verrà scritto sul file *(sequential_output.txt)*. E' stato poi realizzato un programma *(output_correctness.c)* che mette a confronto il contenuto di questi file per verificarne la correttezza.

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

# Problem evaluation and Benchmarks

Per la valutazione delle prestazioni della soluzione proposta sono state create delle versioni dei programmi (sequenziale: sequential_nBody_benchmarking.c - parallela: parallel_nBody_benchmarking.c) leggermente revisionate poichè per effettuare un migliore benchmarking sono state eliminate sia scritture su stdout durante la computazione sia la scrittura dell'output finale su file.

Possiamo ora procedere con la descrizione dei risultati dati dalla  misurazione della scalabilità dell'applicazione. Ci sono due modi di base per misurare la performance parallela di un'applicazione: strong e weak scaling.

### Strong Scaling

In questo tipo di misurazione la taglia del problema (il numero di particelle) resta fissata ma il numero di processori aumenta. 



### Weak Scaling

In questo caso la taglia del problema aumenta con l'aumentare del numero di processori, facendo in modo che il workload sia sempre equamente distribuito tra i processori.





