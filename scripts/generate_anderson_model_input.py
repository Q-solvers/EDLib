#!/usr/bin/env python

import h5py
import numpy as np


def Kanamori_interaction(l, U_int, J_hund):
    norb = 2*l + 1
    U =  np.zeros((norb,norb,norb,norb), dtype=float)
    for i in range(norb):
        U[i][i][i][i] = U_int
        for j in range(norb):
            if(i != j):
                U[i][j][i][j] = U_int - 2* J_hund
                U[i][j][j][i] = J_hund
                U[i][i][j][j] = J_hund
    return U

l = 2

U = np.real(Kanamori_interaction(l, U_int=2.0, J_hund=0.3))
xmu = 4.5
Eps0 = np.array([[0.,0.], [0.,0.], [0.,0.], [0.,0.], [0.,0.] ])
Vk_ =   [ np.array([ [0.5,  0.5] ]),
         np.array([ [0.5,  0.5] ]),
         np.array([ [0.5,  0.5] ]),
         np.array([ [0.5,  0.5] ]),
         np.array([ [0.5,  0.5] ]) ]
Epsk = np.array([ [-1.0,-1.0],
                  [-1.0,-1.0],
                  [-1.0,-1.0],
                  [-1.0,-1.0],
                  [-1.0,-1.0] ])

ml = 2*l + 1
H0 = np.zeros((ml,ml,2))
Vk = np.zeros((ml,) + Epsk.shape)
for i in range(ml):
    H0[i, i,:] = Eps0[i,:]
    Vk[i, i, :] = Vk_[i]

Nk = 1
Ns = Eps0.shape[0] + Epsk.shape[0]

sectors = np.array([[3,3],])

data = h5py.File("input.h5", "w");

beta = data.create_dataset("BETA", shape=(), dtype='f', data=10.0)

hop_g = data.create_group("sectors")
hop_g.create_dataset("values", data=sectors)

bath = data.create_group("Bath")

bath["Epsk/values"] = Epsk
for i in range(ml):
    if(Epsk.shape != Vk[i].shape):
        raise "Incorrect shape for Hybridisation and Epsk"
    Vk_g = bath.create_group("Vk_" + str(i))
    Vk_g.create_dataset("values", data=np.array(Vk[i]), dtype=np.float64)
    t0_g = data.create_group("H0_" + str(i))
    t0_g.create_dataset("values", data=np.array(H0[i]), dtype=np.float64)

UU = np.zeros((2,2) + U.shape)
UU[0,0] = U
UU[0,1] = U
UU[1,0] = U
UU[1,1] = U

int_g = data.create_group("interaction")
int_ds = int_g.create_dataset("values", shape=(2,2,ml,ml,ml,ml,), data=UU)

#int_g = data.create_group("chemical_potential")
int_ds = data.create_dataset("mu", shape=(), data=xmu)

