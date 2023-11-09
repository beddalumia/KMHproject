import numpy as np
import matplotlib.pyplot as plt
import strawberrypy_tmp as sbp
import os
import glob

class cd:
    """Context manager for changing the current working directory"""
    def __init__(self, newPath):
        self.newPath = os.path.expanduser(newPath)

    def __enter__(self):
        self.savedPath = os.getcwd()
        os.chdir(self.newPath)

    def __exit__(self, etype, value, traceback):
        os.chdir(self.savedPath)

for radius in sorted(glob.glob('R=*/')):

    with cd(radius):

        for uloc in sorted(glob.glob('U=*/')):

            with cd(uloc):

                # Tell the user where we are
                print(radius+uloc)

                try:

                    # Build the spin matrix
                    flake = np.loadtxt('flake.txt')
                    Nsite = np.size(flake)
                    np.savetxt("szmatrix.dat", np.diag(np.array([np.array([1, -1]) for _ in range(int(Nsite/2))]).flatten()))

                    # Build finite model from file
                    model = sbp.FiniteModel(mode = 'load',                  # Build from file, not from PythTB or TBmodels
                                spinful = True,                             # Kane-Mele is spinful, needed for counting
                                atoms_uc = 4,                               # Number of states per unit cell (spinful case is 2 * #atoms)
                                uc_vol = 0.8660254037844386,                # Unit cell area (of a Kane-Mele)
                                fnames = ["topoHij.dat",                    # Hamiltonian input file
                                    "spin_xvals.dat",                       # x-positions input file
                                    "spin_yvals.dat",                       # y-positions input file
                                    "szmatrix.dat"                          # Sz matrix input file
                                ])

                    # Evaluate the spin Chern marker
                    spinchern = model.local_spin_chern_marker(macroscopic_average = False, cutoff = 0.9, check_gap = True, histogram = True)

                    # Save to file
                    np.savetxt('Z2marker.txt',spinchern)

                    # Plot the results
                    fig, ax = plt.subplots(1, 1)
                    cbr = ax.scatter(x = model.positions[0], y = model.positions[1], c = spinchern, cmap = "twilight", vmin = 0, vmax = 2)
                    ax.set_title("Topological hamiltonian")
                    plt.colorbar(cbr, label = "Spin Chern marker")
                    ax.set_xlabel("x")
                    ax.set_ylabel("y")
                    ax.set_aspect("equal")
                    plt.tight_layout()
                    fig.savefig("./Z2marker.pdf")
                    plt.close(fig)

                except:
                    print("WARNING: directory skipped, unconverged?")
                else:
                    print("...DONE!")