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

La soluzione proposta ha una complessità quadratica O(n^2) rispetto alla taglia dell'input (il numero di particelle). Esiste un tipo di simulazione, detta di Barnes – Hut che è più efficiente poichè tramite approssimazioni riesce ad eseguire con ordine O(n log n) ma è più complessa da implementare perciò si è preferito considerare la soluzione di complessità quadratica più semplice e concentrarsi sui concetti della programmazione parallela.

# Descrizione della struttura della soluzione
The program must be able to simulate the process for a particular number of steps I. The MASTER process initializes an array of bodies at random and sends it to P-1 processors. Notice that the MASTER can contribute to the computation or not; it is your choice. Each slave simulates the bodies force, only for its bodies, and sends results of its bodies, needed for the next simulation step—the hard part of the problem concerning how to reduce the communication overhead. For the initialization phase, you can consider creating the bodies in a file, and all processors start by reading this file or use a deterministic algorithm to randomize the bodies initialization. The important part is that the execution using one or P processors must generate the same output.

