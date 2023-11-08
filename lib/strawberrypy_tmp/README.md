# StraWBerryPy - Solo marker locali
Ciao! Questa è la mia parte di codice che prima o poi comparirà in strawberrypy.
La struttura cambierà un po' per interfacciare nel modo più semplice possibile la parte a singolo punto di Roberta, ma il funzionamento dovrebbe rimanere simile. In particolare anche per i file di input stiamo anora vedendo che formato chiedere, se i file di wannier90 o file di testo o altro.
Io non sapendo nulla di wannier90 qui ho implementato la versione file txt come ti raccontavo, ma non posso assicurare che rimarrà (se rimarrà diventerà più consistente sul modo in cui le cose vengono lette, qui è fatto un po' a comodità).
Inoltre, la versione che ho qui _dovrebbe_ funzionare correttamente con il formato di file che mi hai mandato per mail (in caso volessi modificare, la lettura dei file avviene in [package.py](strawberrypy_tmp/package.py) dentro la classe ```FiniteModel```, nei metodi ```_load_XXX```).

Ti lascio qui sotto il codice che ho usato per produrre i grafici che ti ho mandato (ho messo i file dentro una cartella "input"):
```python
import numpy as np
import matplotlib.pyplot as plt
import strawberrypy_tmp as sbp

# Build the spin matrix
if True: np.savetxt("./input/szmatrix.dat", np.diag(np.array([np.array([1, -1]) for _ in range(150)]).flatten()))

# Build finite model from file
model = sbp.FiniteModel(mode = 'load',                          # Build from file, not from PythTB or TBmodels
            spinful = True,                                     # Kane-Mele is spinful, needed for counting
            atoms_uc = 4,                                       # Number of states per unit cell (spinful case is 2 * #atoms)
            uc_vol = 0.8660254037844386,                        # Unit cell area (of a Kane-Mele)
            fnames = ["./input/topoHij_noninteracting.dat",     # Hamiltonian input file
                "./input/spin_xvals.dat",                       # x-positions input file
                "./input/spin_yvals.dat",                       # y-positions input file
                "./input/szmatrix.dat"                          # Sz matrix input file
            ])

# Evaluate the spin Chern marker
spinchern = model.local_spin_chern_marker(macroscopic_average = True, cutoff = 0.9, check_gap = True, histogram = True)

# Plot the results
fig, ax = plt.subplots(1, 1)
cbr = ax.scatter(x = model.positions[0], y = model.positions[1], c = spinchern, cmap = "twilight", vmin = 0, vmax = 2)
ax.set_title("Noninteracting hamiltonian")
plt.colorbar(cbr, label = "Spin Chern marker")
ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_aspect("equal")
plt.tight_layout()
fig.savefig("./spinchern.pdf")
```
Ho messo la media macroscopica perchè mi sembra il metodo più comodo per tenere conto del fatto che indici successivi possono non essere vicini, come raccontavi. In questo modo la traccia viene fatta in base alle posizioni dei siti (c'è una funzione di "contrazione" in fondo a [package.py](strawberrypy_tmp/package.py) fa questo).

## Creare l'istanza della classe
La creazione dell'istanza del modello si può fare in due modi, specificati da ```mode```. Se ```mode='tb'``` il la classe si aspetta in input ```tbmodel``` (un modello che arriva da ```PythTB``` o ```TBmodels```), ```nx_sites``` e ```ny_sites``` (il numero di celle unitarie ripetute nelle direzioni _x_ e _y_), e ```spinful``` (```True``` o ```False``` a seconda che il modello abbia o meno spin).
Se invece si vuole leggere il "modello" da file, è necessario specificare ```mode='load'```, ```fnames``` (ovvero il percorso per i file di input, in ordine [matrice hamiltoniana, lista posizioni x, lista posizioni y, matrice spin sz]), ```spinful``` (come nel caso sopra), ```uc_vol``` (area della cella unitaria), ```atoms_uc``` (numero di stati per cella unitaria).

## Local spin Chern marker
Per quanto riguarda il metodo ```local_spin_chern_marker``` questo può prendere diversi argomenti. ```return_projector``` è un bool che se ```True``` ritorna, oltre al marker, una lista con i proiettori [P+, P-] (utile se c'è necessità di fare conti diversi per lo stesso sistema). Se si hanno a dispoizione i proiettori, si possono quindi passare in formato [P+, P-] a ```projector```.
Se il sistema è disordinato o caricato da file, le posizioni dei siti non corrispondono all'indicizzazione interna degli stati, di conseguenza il modo corretto per valutare il marker è fare una media macroscopica (tramite ```macroscopic_average```) che prende in considerazione le posizioni come vengono fornite, e fa quindi una media entro un certo raggio di cutoff (```cutoff```).
Se vuoi controllare se lo spettro di PSzP si chiude o meno ci sono due modi. Uno è il flag ```check_gap``` che controlla priorio il gap tra i due autovalori centrali, e l'altro è l flag ```histogram``` che produce un pdf "pszp_histogram.pdf" della distribuzione degli autovalori di PSzP.
Due argomenti che non ho nominato sono ```direction``` e ```start``` perchè si basano appunto sull'indicizzazione interna, e quindi nel caso di posizioni non ordinate non hanno molto senso.


Credo di aver detto tutto, ma per qualsiasi questione o dubbio scrivimi pure quando vuoi!