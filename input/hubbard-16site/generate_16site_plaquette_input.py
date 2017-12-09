#!/usr/bin/env python

import h5py
import numpy as np

U = np.array([6.0, 6.0, 6.0, 6.0, 6.0, 6.0, 6.0, 6.0, 6.0, 6.0, 6.0, 6.0, 6.0, 6.0, 6.0, 6.0])
xmu = np.array([0.54, 0.54, 0.54, 0.54, 0.54, 0.54, 0.54, 0.54, 0.54, 0.54, 0.54, 0.54, 0.54, 0.54, 0.54, 0.54])

t = -1.0;
tp = 0.3;

sectors = np.array([[6,6]])


Nrows = 4;
Ncols = 4;

tm = np.zeros((Nrows*Ncols, Nrows*Ncols));

for row in range(0, Nrows):
  for col in range(0, Ncols):
    i = row * Ncols + col;
    # Horizontal hopping with period.
    j = row * Ncols + (col + 1)%Ncols;
    tm[i, j] = -t;
    tm[j, i] = -t;
    # Vertical hopping with period.
    j = ((row + 1)%Nrows) * Ncols + col;
    tm[i, j] = -t;
    tm[j, i] = -t;
    # Diagonal (\) hopping with period.
    j = ((row + 1)%Nrows) * Ncols + (col + 1)%Ncols;
    tm[i, j] = -tp;
    tm[j, i] = -tp;
    # Diagonal (/) hopping with period.
    j = ((row + 1)%Nrows) * Ncols + (col - 1)%Ncols;
    tm[i, j] = -tp;
    tm[j, i] = -tp;

print "bandwidth: %f" % (np.amax(np.linalg.eig(tm)[0]) - np.amin(np.linalg.eig(tm)[0]));


Ns = len(U)

data = h5py.File("input.h5", "w");

hop_g = data.create_group("sectors")
hop_g.create_dataset("values", data=sectors)


hop_g = data.create_group("hopping")
hop_g.create_dataset("values", data=tm)

int_g = data.create_group("interaction")
int_ds = int_g.create_dataset("values", shape=(Ns,), data=U)

int_g = data.create_group("chemical_potential")
int_ds = int_g.create_dataset("values", shape=(Ns,), data=xmu)
