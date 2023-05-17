#!/usr/bin/env python

# Two dimensional tight-binding 2D Kane-Mele model
# C.L. Kane and E.J. Mele, PRL 95, 146802 (2005) Eq. (1)

# Copyright under GNU General Public License 2010, 2012, 2016
# by Sinisa Coh and David Vanderbilt (see gpl-pytb.txt)

from __future__ import print_function
from pythtb import * # import TB model class
import numpy as np
import matplotlib.pyplot as plt

def get_kane_mele(t1,t2,tr,mh):
  "Return a Kane-Mele model in the normal or topological phase."

  # define lattice vectors
  lat=[[1.0,0.0],[0.5,np.sqrt(3.0)/2.0]]
  # define coordinates of orbitals
  orb=[[1./3.,1./3.],[2./3.,2./3.]]
  
  # make two dimensional tight-binding Kane-Mele model
  ret_model=tb_model(2,2,lat,orb,nspin=2)
  
  # set parameters of the model
  thop=t1
  esite=mh
  rashba=tr
  spin_orb=t2  
  
  # set on-site energies
  ret_model.set_onsite([esite,(-1.0)*esite])
  
  # pauli matrices
  sigma_x=np.array([0.,1.,0.,0])
  sigma_y=np.array([0.,0.,1.,0])
  sigma_z=np.array([0.,0.,0.,1])
  
  ###############################################################################
  # set hoppings (one for each connected pair of orbitals)                      #
  # > ret_model.set_hop(amplitude, i, j, [lattice vector to cell containing j]) #
  ###############################################################################
  
  # spin-independent first-neighbor hoppings
  ret_model.set_hop(thop, 0, 1, [ 0, 0])
  ret_model.set_hop(thop, 0, 1, [ 0,-1])
  ret_model.set_hop(thop, 0, 1, [-1, 0])
  
  # second-neighbour spin-orbit hoppings (s_z)
  ret_model.set_hop(-1.j*spin_orb*sigma_z, 0, 0, [ 0, 1])
  ret_model.set_hop( 1.j*spin_orb*sigma_z, 0, 0, [ 1, 0])
  ret_model.set_hop(-1.j*spin_orb*sigma_z, 0, 0, [ 1,-1])
  ret_model.set_hop( 1.j*spin_orb*sigma_z, 1, 1, [ 0, 1])
  ret_model.set_hop(-1.j*spin_orb*sigma_z, 1, 1, [ 1, 0])
  ret_model.set_hop( 1.j*spin_orb*sigma_z, 1, 1, [ 1,-1])
  
  # Rashba first-neighbor hoppings: (s_x)(dy)-(s_y)(d_x)
  r3h =np.sqrt(3.0)/2.0
  # bond unit vectors are (r3h,half) then (0,-1) then (-r3h,half)
  ret_model.set_hop(1.j*rashba*( 0.5*sigma_x-r3h*sigma_y), 0, 1, [ 0, 0], mode="add")
  ret_model.set_hop(1.j*rashba*(-1.0*sigma_x            ), 0, 1, [ 0,-1], mode="add")
  ret_model.set_hop(1.j*rashba*( 0.5*sigma_x+r3h*sigma_y), 0, 1, [-1, 0], mode="add")

  return ret_model

# now solve the model and find Wannier centers for both topological
# and normal phase of the model
for t1 in [1]:
	for t2 in [0.01,0.1,0.3,0.55]:
		for tr in [0]:
			for mh in [0,1.6]:
				
				# get the tight-binding model
				my_model=get_kane_mele(t1,t2,tr,mh)

				# list of nodes (high-symmetry points) that will be connected
				path=[[0.,0.],[2./3.,1./3.],[.5,.5],[1./3.,2./3.], [0.,0.]]
				# labels of the nodes
				label=(r'$\Gamma $',r'$K$', r'$M$', r'$K^\prime$', r'$\Gamma $')
				(k_vec,k_dist,k_node)=my_model.k_path(path,101,report=False)
				
				# initialize figure with subplots
				fig, (ax1, ax2) = plt.subplots(1,2,figsize=(6.5,2.8))
				
				# solve for eigenenergies of hamiltonian on
				# the set of k-points from above
				evals=my_model.solve_all(k_vec)
				# plot bands
				ax1.plot(k_dist,evals[0],color='k')
				ax1.plot(k_dist,evals[1],color='k')
				ax1.plot(k_dist,evals[2],color='k')
				ax1.plot(k_dist,evals[3],color='k')
				ax1.set_title("$\lambda_{SO}=$"+str(t2)+"t and $M_h=$"+str(mh))
				ax1.set_xticks(k_node)
				ax1.set_xticklabels(label)
				ax1.set_xlim(k_node[0],k_node[-1])
				for n in range(len(k_node)):
					ax1.axvline(x=k_node[n],linewidth=0.5, color='k')
				ax1.set_xlabel("k-space")
				ax1.set_ylabel("Energy")
				ax1.set_ylim(-3,3)
				
				#calculate my-array
				my_array=wf_array(my_model,[41,41])
				
				# solve model on a regular grid, and put origin of
				# Brillouin zone at [-1/2,-1/2]  point
				my_array.solve_on_grid([-0.5,-0.5])
				
				# compute chern numbers by integrating berry curvature
				print("Chern number for filled spin-bands")
				print("  C = %5.2f" % ((my_array.berry_flux([0])+my_array.berry_flux([1]))/(2.*np.pi)))
				print("")

				# calculate Berry phases around the BZ in the k_x direction
				# (which can be interpreted as the 1D hybrid Wannier centers
				# in the x direction) and plot results as a function of k_y
				#
				# Following the ideas in
				#   A.A. Soluyanov and D. Vanderbilt, PRB 83, 235401 (2011)
				#   R. Yu, X.L. Qi, A. Bernevig, Z. Fang and X. Dai, PRB 84, 075119 (2011)
				# the connectivity of these curves determines the Z2 index
				#
				wan_cent = my_array.berry_phase([0,1],dir=1,contin=False,berry_evals=True)
				wan_cent/=(2.0*np.pi)
				
				nky=wan_cent.shape[0]
				ky=np.linspace(0.,1.,nky)
				# draw shifted Wannier center positions
				for shift in range(-2,3):
					ax2.plot(ky,wan_cent[:,0]+float(shift),"k.")
					ax2.plot(ky,wan_cent[:,1]+float(shift),"k.")
				ax2.set_ylim(-1.0,1.0)
				ax2.set_ylabel('Wannier center along x')
				ax2.set_xlabel(r'$k_y$')
				ax2.set_xticks([0.0,0.5,1.0])
				ax2.set_xlim(0.0,1.0)
				ax2.set_xticklabels([r"$0$",r"$\pi$", r"$2\pi$"])
				ax2.axvline(x=.5,linewidth=0.5, color='k')
				ax2.set_title("1D Wannier centers")

				fig.tight_layout()
				fig.savefig("kane_mele_so"+str(t2)+"_mh"+str(mh)+".svg",format='svg')

print('Done.\n')

