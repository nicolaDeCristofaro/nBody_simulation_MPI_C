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

# Descrizione della struttura della soluzione
The program must be able to simulate the process for a particular number of steps I. The MASTER process initializes an array of bodies at random and sends it to P-1 processors. Notice that the MASTER can contribute to the computation or not; it is your choice. Each slave simulates the bodies force, only for its bodies, and sends results of its bodies, needed for the next simulation step—the hard part of the problem concerning how to reduce the communication overhead. For the initialization phase, you can consider creating the bodies in a file, and all processors start by reading this file or use a deterministic algorithm to randomize the bodies initialization.  **DA COMPLETARE**

# Correttezza
Per la correttezza di un programma in parallelo è necessario che l'esecuzione con P processori o con un solo processore con lo stesso input produca lo stesso output. Per verificare questa condizione è stato fatto in modo che l'output dell'esecuzione parallela con più processori sia stato scritto su file (nome file) così come l'output della versione parallela (nome file). E' stato poi realizzato un programma (nome file programma) che mette a confronto il contenuto di questi file per verificarne la correttezza.

Dopo aver eseguito la versione parallela e quella sequenziale -> comando per eseduire correttezza

Da sottolineare che per quanto riguarda i confronti tra gli attributi delle particelle per stabilirne la correttezza è stata utilizzata una funzione per confrontare i valori float. Questo perchè la matematica in virgola mobile non è esatta. Valori semplici come 0,2 non possono essere rappresentati con precisione usando numeri binari in virgola mobile e questa precisione limitata dei numeri in virgola mobile significa che lievi variazioni nell'ordine delle operazioni possono cambiare il risultato. Quindi, dato che questi valori vengono memorizzati con precisioni diverse, i risultati possono differire. Di conseguenza se si eseguono vari calcoli e quindi si confrontano i risultati con un valore atteso, è altamente improbabile che si ottenga esattamente il risultato desiderato. Per questo motivo con la funzione realizzata possiamo esprimere il concetto che due valori abbastanza vicini tra loro possono essere considerati uguali.

In pratica, due valori di tipo float sono considerati uguali se la loro differenza rientra in un certo limite o valore epsilon.

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

