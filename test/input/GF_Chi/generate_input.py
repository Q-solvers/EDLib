#!/usr/bin/env python

import h5py
import numpy as np
np.set_printoptions(linewidth=999)

U = 3.0
xmu = 1.0

t = -1.0
tp = 0.3

gforb = np.array([
 [0,0]
])
chiorb = np.array([
 [0,0]
])

Nrows = 2
Ncols = 2

Um = np.zeros(Nrows * Ncols)
xmum = np.zeros(Nrows * Ncols)

for row in range(0, Nrows):
  for col in range(0, Ncols):
   ii = row * Ncols + col
   Um[ii] = U
   xmum[ii] = xmu

tm = np.zeros((Nrows*Ncols, Nrows*Ncols))

for row in range(0, Nrows):
  for col in range(0, Ncols):
    i = row * Ncols + col
    # Vertical hopping.
    j = ((row + 1)%Nrows) * Ncols + col
    tm[i, j] = -t
    tm[j, i] = -t
    # Horizontal hopping.
    j = row * Ncols + (col + 1)%Ncols
    tm[i, j] = -t
    tm[j, i] = -t
    # Diagonal (\) hopping.
    j = ((row + 1)%Nrows) * Ncols + (col + 1)%Ncols
    tm[i, j] = -tp
    tm[j, i] = -tp
    # Diagonal (/) hopping.
    j = ((row + 1)%Nrows) * Ncols + (col - 1)%Ncols
    tm[i, j] = -tp
    tm[j, i] = -tp

print "U:", Um
print "xmu:", xmum
print "t:"
print tm
print "bandwidth: %f" % (np.amax(np.linalg.eig(tm)[0]) - np.amin(np.linalg.eig(tm)[0]))


Ns = len(Um)

data = h5py.File("input.h5", "w")


hop_g = data.create_group("hopping")
hop_g.create_dataset("values", data=tm)

int_g = data.create_group("interaction")
int_ds = int_g.create_dataset("values", shape=(Ns,), data=Um)

chem_g = data.create_group("chemical_potential")
chem_ds = chem_g.create_dataset("values", shape=(Ns,), data=xmum)

gf_g = data.create_group("GreensFunction_orbitals")
gf_g.create_dataset("values", data=gforb)

chi_g = data.create_group("ChiLoc_orbitals")
chi_g.create_dataset("values", data=chiorb)
