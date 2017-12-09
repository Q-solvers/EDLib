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
Vk =   [ np.array([ [0.5,  0.5] ]),
         np.array([ [0.5,  0.5] ]),
         np.array([ [0.5,  0.5] ]),
         np.array([ [0.5,  0.5] ]),
         np.array([ [0.5,  0.5] ]) ]
Epsk = [ np.array([ [-1.0,-1.0] ]),
         np.array([ [-1.0,-1.0] ]),
         np.array([ [-1.0,-1.0] ]),
         np.array([ [-1.0,-1.0] ]),
         np.array([ [-1.0,-1.0] ]) ]

tk = [ np.array([[0, 0], [0, 0], [0, 0], [0, 0], [0, 0]]),
       np.array([[0, 0], [0, 0], [0, 0], [0, 0], [0, 0]]),
       np.array([[0, 0], [0, 0], [0, 0], [0, 0], [0, 0]]),
       np.array([[0, 0], [0, 0], [0, 0], [0, 0], [0, 0]]),
       np.array([[0, 0], [0, 0], [0, 0], [0, 0], [0, 0]]) ]

ml = 2*l + 1
Nk = 1
Ns = len(Eps0) + len(Epsk)

sectors = np.array([[3,3],])

data = h5py.File("input.h5", "w");

beta = data.create_dataset("BETA", shape=(), dtype='f', data=10.0)

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
    t0_g.create_dataset("values", data=np.array(tk[i]), dtype=np.float)

hop_g = data.create_group("Eps0")
hop_g.create_dataset("values", data=Eps0)

int_g = data.create_group("interaction")
int_ds = int_g.create_dataset("values", shape=(ml,ml,ml,ml,), data=U)

#int_g = data.create_group("chemical_potential")
int_ds = data.create_dataset("mu", shape=(), data=xmu)

