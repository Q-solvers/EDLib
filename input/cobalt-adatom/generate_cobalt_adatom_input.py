#!/usr/bin/env python

import h5py
import numpy as np

# We use pytriqs package for coulomb matrix calculation
import pytriqs.operators.util as U_matrix


# Coulomb matrix
U = np.real(U_matrix.U_matrix(2, U_int=6.6, J_hund=0.9, basis='cubic'))
# Chemical potential
xmu = 44.44

# Impurity energies
Eps0 = np.array([[-0.70,-0.70], [-0.8,-0.8], [-0.50,-0.50], [-0.8,-0.8], [-0.70,-0.70] ])
# Hybridization between impurity and bath
Vk =   [ np.array([ [0.56434,0.56434],[0.68392,0.68392],[0.29519,0.29519]]),
         np.array([ [0.81892,0.81892],[0.99136,0.99136]]),
         np.array([ [0.77347,0.77347],[0.79785,0.79785]]),
         np.array([ [0.81892,0.81892],[0.99136,0.99136]]),
         np.array([ [0.56434,0.56434],[0.68392,0.68392],[0.29519,0.29519]]) ]
# Bath energies
Epsk = [ np.array([ [-2.37325, -2.37325],[-0.87328,-0.87328], [2.01265,2.01265]]),
         np.array([ [-3.15496,-3.15496],[-1.69066,-1.69066] ]),
         np.array([ [-5.59842,-5.59842],[-2.95325,-2.95325]]),
         np.array([ [-3.15496,-3.15496],[-1.69066,-1.69066]]),
         np.array([ [-2.37325, -2.37325],[-0.87328,-0.87328], [2.01265,2.01265]] ) ]

# Total number of fermionic orbitals
Ns = 17
# Number of impurity orbitals
ml = 5

# Interorbital hoppings
t0 = [ np.array([[0, 0]]*ml),
       np.array([[0, 0]]*ml),
       np.array([[0, 0]]*ml),
       np.array([[0, 0]]*ml),
       np.array([[0, 0]]*ml) ]

# (n_up, n_down) symmerty sectors
sectors = np.array([[12,15],[15,12],[14,13],[13,14],[13,15],[14,14],[15,13],])

# open data file
data = h5py.File("input.h5", "w");
# Inverse temperature
#beta = data.create_dataset("BETA", shape=(), dtype='f', data=10.0)

hop_g = data.create_group("sectors")
hop_g.create_dataset("values", data=sectors)

bath = data.create_group("Bath")

for i in range(ml):
  if(Epsk[i].shape != Vk[i].shape):
    raise "Incorrect shape for Hybridisation and Epsk"
  Epsk_g = bath.create_group("Epsk_" + str(i))
  Epsk_g.create_dataset("values", shape=(len(Epsk[i]),2,), data=Epsk[i], dtype=np.float)
  Vk_g = bath.create_group("Vk_" + str(i))
  Vk_g.create_dataset("values", data=np.array(Vk[i]), dtype=np.float)
  t0_g = data.create_group("t0_" + str(i))
  t0_g.create_dataset("values", data=np.array(t0[i]), dtype=np.float)

hop_g = data.create_group("Eps0")
hop_g.create_dataset("values", data=Eps0)

int_g = data.create_group("interaction")
int_ds = int_g.create_dataset("values", shape=(ml,ml,ml,ml,), data=U)

int_ds = data.create_dataset("mu", shape=(), data=xmu)

