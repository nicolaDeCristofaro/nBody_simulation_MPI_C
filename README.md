| N-Body Simulation | Nicola De Cristofaro | Data di consegna |
| --- | --- | --- |

# Descrizione del problema
In fisica, il problema n-body consiste nel predire i  singoli movimenti di un gruppo di oggetti celesti che interagiscono tra loro in modo gravitazionale. La risoluzione di questo problema è stata motivata dal desiderio di comprendere i movimenti del Sole, della Luna, dei pianeti e delle stelle visibili.

<div style="text-align:center;display: block;"><img alt="nBodySimulation" src="https://github.com/nicolaDeCristofaro/nBody_simulation_MPI_C/blob/master/images/NaturalDiligentIguanodon-size_restricted.gif" /></div>

# Descrizione della soluzione
Per trovare la soluzione a questo problema è possibile simulare il comportamento delle particelle, ognuna avente una massa, una posizione e una velocità iniziale. La simulazione consentirà di calcolare la posizione e la velocità di ciascuna particella in una sequenza di intervalli di tempo specifici, o semplicemente la posizione e la velocità di ciascuna particella al termine di intervallo di tempo determinato.

La soluzione proposta ha una complessità quadratica O(n^2) rispetto alla taglia dell'input (il numero di particelle).Esiste un tipo di simulazione, detta di Barnes – Hut che è più efficiente poichè tramite approssimazioni riesce ad eseguire con ordine O(n log n) ma è più complessa da implementare perciò si è preferito considerare la soluzione di complessità quadratica più semplice e concentrarsi sui concetti della programmazione parallela.
