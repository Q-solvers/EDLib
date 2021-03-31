#!/usr/bin/env python

import h5py
import numpy as np
np.set_printoptions(linewidth=999)

H = 0
checkerboard_afm = 0

U = 6.0
xmu = 0.0

t = -1.0
tp = 0.3

nonperiodic = 0

alpha = 1.0

at = alpha * t
atp = alpha * tp

sectors = np.array([
 [7,7],
 [8,8]
])

gforb = np.array([
])
chiorb = np.array([
])


Nrows = 4
Ncols = 4

Um = np.zeros(Nrows * Ncols)
xmum = np.zeros(Nrows * Ncols)

for row in range(0, Nrows):
  for col in range(0, Ncols):
   ii = row * Ncols + col
   Um[ii] = U
   xmum[ii] = xmu

tm = np.zeros((Nrows*Ncols, Nrows*Ncols))

# Vertical hopping.
for row in range(0, Nrows - nonperiodic):
  for col in range(0, Ncols):
    i = row * Ncols + col
    j = ((row + 1)%Nrows) * Ncols + col
    if ((row % 2) == 0):
     # t
     tm[i, j] += -t
     tm[j, i] += -t
    else:
     # alpha t
     tm[i, j] += -at
     tm[j, i] += -at

# Horizontal hopping.
for row in range(0, Nrows):
  for col in range(0, Ncols - nonperiodic):
    i = row * Ncols + col
    j = row * Ncols + (col + 1)%Ncols
    if ((col % 2) == 0):
     # t
     tm[i, j] += -t
     tm[j, i] += -t
    else:
     # alpha t
     tm[i, j] += -at
     tm[j, i] += -at

# Diagonal (\) hopping.
for row in range(0, Nrows - nonperiodic):
  for col in range(0, Ncols - nonperiodic):
    i = row * Ncols + col
    j = ((row + 1)%Nrows) * Ncols + (col + 1)%Ncols
    if (((col % 2)  ==  0) & ((row % 2)  ==  0)):
     # tp
     tm[i, j] += -tp
     tm[j, i] += -tp
    else:
     # alpha tprime
     tm[i, j] += -atp
     tm[j, i] += -atp

# Diagonal (/) hopping.
for row in range(0, Nrows - nonperiodic):
  for col in range(nonperiodic, Ncols):
    i = row * Ncols + col
    j = ((row + 1)%Nrows) * Ncols + (col - 1)%Ncols
    if (((col % 2) != 0) & ((row % 2) == 0)):
     # tp
     tm[i, j] += -tp
     tm[j, i] += -tp
    else:
     # alpha tprime
     tm[i, j] += -atp
     tm[j, i] += -atp

Hm = np.zeros(Nrows * Ncols)
for row in range(0, Nrows - nonperiodic):
  for col in range(nonperiodic, Ncols):
   site = row * Ncols + col
   if(checkerboard_afm):
    Hm[site] = (1.0 - ((row + col) % 2) * 2.0) * H
   else:
    Hm[site] = H

print("U:", Um)
print("xmu:", xmum)
print("H:", Hm)
print("t:")
print(tm)
print("bandwidth: %f" % (np.amax(np.linalg.eig(tm)[0]) - np.amin(np.linalg.eig(tm)[0])))


Ns = len(Um)

data = h5py.File("input.h5", "w")

hop_g = data.create_group("sectors")
hop_g.create_dataset("values", data=sectors)


hop_g = data.create_group("hopping")
hop_g.create_dataset("values", data=tm)

int_g = data.create_group("interaction")
int_ds = int_g.create_dataset("values", shape=(Ns,), data=Um)

mag_g = data.create_group("magnetic_field")
mag_ds = mag_g.create_dataset("values", shape=(Ns,), data=Hm)

chem_g = data.create_group("chemical_potential")
chem_ds = chem_g.create_dataset("values", shape=(Ns,), data=xmum)

gf_g = data.create_group("GreensFunction_orbitals")
gf_g.create_dataset("values", data=gforb)

chi_g = data.create_group("ChiLoc_orbitals")
chi_g.create_dataset("values", data=chiorb)
