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

import logging

import z2pack
import numpy as np

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

def Hk_kanemele(k, Mh, t1, t2, phi):
    # LATTICE VECTORS FOR A HONEYCOMB
    # e₁ = a₀ [ sqrt3/2 , 1/2 ] = 3/2a[1, 1/sqrt3]
    # e₂ = a₀ [ sqrt3/2 ,-1/2 ] = 3/2a[1,-1/sqrt3]
    e1 = 3/2. * np.array([1, 1/np.sqrt(3)])
    e2 = 3/2. * np.array([1,-1/np.sqrt(3)])
    # 
    kdote1 = np.dot(k,e1)
    kdote2 = np.dot(k,e2)
    #
    h0 = 2*t2*np.cos(phi)*(np.cos(kdote1) + np.cos(kdote2) + np.cos(kdote1-kdote2))
    hx = t1*(np.cos(kdote1) + np.cos(kdote2) + 1)
    hy = t1*(np.sin(kdote1) + np.sin(kdote2))
    hz = 2*t2*np.sin(phi)*(np.sin(kdote1) - np.sin(kdote2) - np.sin(kdote1-kdote2))
    #
    HkUP = h0*pauli_0 + hx*pauli_x + hy*pauli_y + hz*pauli_z + Mh*pauli_z
    HkDW = h0*pauli_0 + hx*pauli_x + hy*pauli_y - hz*pauli_z + Mh*pauli_z
    #
    Hk = np.block([
        [HkUP           ,   np.zeros((2,2))],
        [np.zeros((2,2)),   HkDW]
         ])
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
        lambda k: Hk_kanemele(k, m, t1, t2, phi), dim=2, bands=2
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
        system=system_UP, surface=lambda s, t: [s,t], min_neighbour_dist=1e-5
    )
    result_DW = z2pack.surface.run(
        system=system_DW, surface=lambda s, t: [s,t], min_neighbour_dist=1e-5
    )

    Z2_I = z2pack.invariant.z2(result)
    print(" Z2 [direct evaluation] = %5.2f" % Z2_I)
    C_UP = z2pack.invariant.chern(result_UP)
    C_DW = z2pack.invariant.chern(result_DW)
    print(" Z2 [ (Cup - Cdw) / 2 ] = %5.2f" % ((C_UP - C_DW) / 2.0))
    print("")

    return (C_UP - C_DW) / 2.0


if __name__ == "__main__":
    for t1 in [1]:
        for t2 in [0.01,0.1,0.3,0.55,1]:
            for mh in [0,1.6]:
                get_invariant(mh, t1, t2, 0.5 * np.pi)