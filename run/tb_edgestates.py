#!/usr/bin/env python

# Haldane model from Phys. Rev. Lett. 61, 2015 (1988)
# Solves model and draws one of its edge states.

# Copyright under GNU General Public License 2010, 2012, 2016
# by Sinisa Coh and David Vanderbilt (see gpl-pythtb.txt)

from __future__ import print_function
from pythtb import * # import TB model class
import numpy as np

# define lattice vectors
lat=[[3.0/2.0,-np.sqrt(3)/2.0],[3.0/2.0,np.sqrt(3.0)/2.0]]

# define coordinates of orbitals
orb=[[1./3.,1./3.],[2./3.,2./3.]]

# make two dimensional tight-binding Haldane model
my_model=tb_model(2,2,lat,orb)

# set model parameters
delta=0.0
t=-1.0
t2=0.01*np.exp((1.j)*np.pi/2.)
t2c=t2.conjugate()

# set on-site energies
my_model.set_onsite([-delta,delta])
# set hoppings (one for each connected pair of orbitals)
# (amplitude, i, j, [lattice vector to cell containing j])
my_model.set_hop(t, 0, 1, [ 0, 0])
my_model.set_hop(t, 1, 0, [ 1, 0])
my_model.set_hop(t, 1, 0, [ 0, 1])
# add second neighbour complex hoppings
my_model.set_hop(t2 , 0, 0, [ 1, 0])
my_model.set_hop(t2 , 1, 1, [ 1,-1])
my_model.set_hop(t2 , 1, 1, [ 0, 1])
my_model.set_hop(t2c, 1, 1, [ 1, 0])
my_model.set_hop(t2c, 0, 0, [ 1,-1])
my_model.set_hop(t2c, 0, 0, [ 0, 1])

# print tight-binding model details
my_model.display()

# cutout finite model first along direction x with no PBC
tmp_model=my_model.cut_piece(10,0,glue_edgs=False)
# cutout also along y direction with no PBC
fin_model=tmp_model.cut_piece(10,1,glue_edgs=False)

# solve finite models
(evals,evecs)=fin_model.solve_all(eig_vectors=True)
with open('levels.txt', 'w') as f:
   for i in range(0,len(evals)):
      print(evals[i],file=f)

# pick index of state in the middle of the gap
ed=fin_model.get_num_orbitals()//2

# pick index of single-particle ground-state
gs=0

# draw one of the edge-states 
(fig,ax)=fin_model.visualize(0,1,eig_dr=evecs[ed,:],draw_hoppings=False)
ax.set_title("Edge state for finite model without periodic direction")
ax.set_xlabel("x coordinate")
ax.set_ylabel("y coordinate")
fig.tight_layout()
fig.savefig("edge_state1.pdf")

# draw the other edge-state
(fig,ax)=fin_model.visualize(0,1,eig_dr=evecs[ed+1,:],draw_hoppings=False)
ax.set_title("Edge state for finite model without periodic direction")
ax.set_xlabel("x coordinate")
ax.set_ylabel("y coordinate")
fig.tight_layout()
fig.savefig("edge_state2.pdf")

# draw the single-particle ground-state
(fig,ax)=fin_model.visualize(0,1,eig_dr=evecs[gs,:],draw_hoppings=False)
ax.set_title("Ground state for finite model without periodic direction")
ax.set_xlabel("x coordinate")
ax.set_ylabel("y coordinate")
fig.tight_layout()
fig.savefig("ground_state.pdf")

print('Done.\n')