#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Two-dimensional Kane Mele (Hubbard) model, defined as an
# explicit Hamiltonian within the Z2Pack formalism. The Z2
# invariant is then computed with the library. The Hubbard
# interaction is accounted for renormalizing the H(k) with
# the so called topological self-energy Σ(iω=0), proven to
# give the same exact topological properties of the fully
# interacting (local) problem Σ(iω) as far as there is an
# adiabatic connection to the interacting model. The usual
# mean-field approach would instead renormalize by Σ(ω->∞)
# giving a poorer description of the interaction/topology
# interplay, given that Hartree-Fock theory optimizes the 
# energy of a single particle solution, without any claim
# on preserving the topological nature of the interacting
# system.
#
# See also https://doi.org/10.1088/0953-8984/25/15/155601

import logging

import z2pack
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
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

logging.getLogger('z2pack').setLevel(logging.WARNING)

# defining pauli matrices
pauli_0 = np.identity(2, dtype=complex)
pauli_x = np.array([[0, 1], [1, 0]], dtype=complex)
pauli_y = np.array([[0, -1j], [1j, 0]], dtype=complex)
pauli_z = np.array([[1, 0], [0, -1]], dtype=complex)


def Hk_haldane(k, m, t1, t2, phi):
    kx, ky = k
    k_a = 2 * np.pi / 3. * np.array([-kx - ky, 2. * kx - ky, -kx + 2. * ky])
    k_b = 2 * np.pi * np.array([kx, -kx + ky, ky])
    H = 2 * t2 * np.cos(phi) * sum([np.cos(-val) for val in k_b]) * pauli_0
    H += t1 * sum([np.cos(-val) for val in k_a]) * pauli_x
    H += t1 * sum([np.sin(-val) for val in k_a]) * pauli_y
    H += m * pauli_z
    H -= 2 * t2 * np.sin(phi) * sum([np.sin(-val) for val in k_b]) * pauli_z
    return H

def Hk_dmft(k, Mh, t1, t2, phi):
    # LATTICE VECTORS FOR A HONEYCOMB
    # e₁ = a₀ [ sqrt3/2 , 1/2 ] = 3/2a[1, 1/sqrt3]
    # e₂ = a₀ [ sqrt3/2 ,-1/2 ] = 3/2a[1,-1/sqrt3]
    e1 = 3/2. * np.array([1, 1/np.sqrt(3)])
    e2 = 3/2. * np.array([1,-1/np.sqrt(3)])
    # 
    kdote1 = np.dot(k,e1) * np.pi * 2
    kdote2 = np.dot(k,e2) * np.pi * 2
    #
    h0 = 2*t2*np.cos(phi)*(np.cos(kdote1) + np.cos(kdote2) + np.cos(kdote1-kdote2))
    hx = t1*(np.cos(kdote1) + np.cos(kdote2) + 1)
    hy = t1*(np.sin(kdote1) + np.sin(kdote2))
    hz = 2*t2*np.sin(phi)*(np.sin(kdote1) - np.sin(kdote2) - np.sin(kdote1-kdote2))
    #
    Hk = h0*pauli_0 + hx*pauli_x + hy*pauli_y + hz*pauli_z + Mh*pauli_z
    return Hk

def Hk_kanemele(k, Mh, t1, t2, phi):
    # NONINTERACTING MODEL
    HkUP = Hk_haldane(k, Mh, t1, t2,  phi)
    HkDW = Hk_haldane(k, Mh, t1, t2, -phi)
    HkUD = np.zeros((2,2))
    HkDU = np.zeros((2,2))
    # TOPOLOGICAL EFFECTIVE INTERACTING MODEL
    # retrieve self-energy corrections
    s11 = pd.read_csv ("impSigma_l11_s11_iw.ed",names=["iw","im","re"],sep='\s+',nrows=1)
    s12 = pd.read_csv ("impSigma_l11_s12_iw.ed",names=["iw","im","re"],sep='\s+',nrows=1)
    s21 = pd.read_csv ("impSigma_l11_s21_iw.ed",names=["iw","im","re"],sep='\s+',nrows=1)
    s22 = pd.read_csv ("impSigma_l11_s22_iw.ed",names=["iw","im","re"],sep='\s+',nrows=1)
    # build spin-up self-energy matrix
    sAup = s11["re"].to_numpy()+1j*s11["im"].to_numpy()
    S_UP = np.block([ 
        [sAup,    0],
        [0   , sAup]  # Σ(iω)_{B,up,up} = Σ(iω)_{A,up,up}
    ]) 
    # build spin-dw self-energy matrix
    sAdw = s22["re"].to_numpy()+1j*s22["im"].to_numpy()
    S_DW = np.block([ 
        [sAdw,    0],
        [0   , sAdw]  # Σ(iω)_{B,dw,dw} = Σ(iω)_{A,dw,dw}
    ]) 
    # build cross-spin self-energy matrices
    sAud = s12["re"].to_numpy()+1j*s12["im"].to_numpy()
    S_UD = np.block([ 
        [sAud,     0],
        [0   , -sAud]  # Σ(iω)_{B,up,dw} = -Σ(iω)_{A,up,dw}
    ]) 
    sAdu = s21["re"].to_numpy()+1j*s21["im"].to_numpy()
    S_DU = np.block([ 
        [sAdu,     0],
        [0   , -sAdu]  # Σ(iω)_{B,dw,up} = -Σ(iω)_{A,dw,up}
    ]) 
    # renormalize the hamiltonian
    HkUP = HkUP + S_UP
    HkDW = HkDW + S_DW
    HkUD = HkUD + S_UD
    HkDU = HkDU + S_DU
    # put it all together
    Hk = np.block([
        [HkUP, HkUD],
        [HkDU, HkDW]
    ])
    #
    return Hk

def get_invariant(m, t1, t2, phi):

    print(" ---------------")
    print(" KANE MELE model")
    print(" ---------------")
    print(" t1 =  1.00")
    print(" t2 = %5.2f" % t2)
    print(" Mh = %5.2f" % m )
    print(" φ  = %4.1fπ"% (phi/np.pi))
    print(" ---------------")

    system = z2pack.hm.System(
        lambda k: Hk_kanemele(k, m, t1, t2, phi), dim=2, bands=2, hermitian_tol=None
    )
    system_UP = z2pack.hm.System(
        lambda k: Hk_haldane(k, m, t1, t2,  phi), dim=2, bands=1
    )
    system_DW = z2pack.hm.System(
        lambda k: Hk_haldane(k, m, t1, t2, -phi), dim=2, bands=1
    )

    result = z2pack.surface.run(
        system=system, surface=lambda s, t: [s/2,t], min_neighbour_dist=1e-5
    )
    result_UP = z2pack.surface.run(
        system=system_UP, surface=lambda s, t: [t,s], min_neighbour_dist=1e-5
    )
    result_DW = z2pack.surface.run(
        system=system_DW, surface=lambda s, t: [t,s], min_neighbour_dist=1e-5
    )

    z2pack.plot.wcc(result)
    plt.savefig("wcc_plot_so"+str(t2)+"_mh"+str(mh)+".svg", bbox_inches='tight')
    plt.close
    z2pack.plot.chern(result_UP)
    plt.savefig("cup_plot_so"+str(t2)+"_mh"+str(mh)+".png", bbox_inches='tight')
    plt.close
    Z2_I = z2pack.invariant.z2(result,check_kramers_pairs=False) # this check is driving me crazy
    print(" Z2 [topological hamiltonian] = %5.2f" % Z2_I)
    C_UP = z2pack.invariant.chern(result_UP)
    C_DW = z2pack.invariant.chern(result_DW)
    print(" Z2 [non-interacting model  ] = %5.2f" % ((C_UP - C_DW) / 2.0))
    print("")

    return Z2_I


if __name__ == "__main__":
    for t1 in [1]:
        for t2 in [0.3]:
            for mh in [0]:
                for uloc in glob.glob('U=*/'):
                    with cd(uloc):
                        print("%s" % uloc)
                        Z2 = get_invariant(mh, t1, t2, 0.5 * np.pi)
                        with open("Z2_inv.dat", "w") as f:
                            print("%f" % Z2, file=f)
