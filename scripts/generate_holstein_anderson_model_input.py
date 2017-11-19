#!/usr/bin/env python

import h5py
import numpy as np


avg = 1.0
U = 3.0
V = 0.0
xmu = 1.5
Eps0 = np.array([[0.,0.]])

# inter-orbital hoppings
tk =   np.array([ [ 0 ] ] )

Vk =   np.array([ [ [0.5,0.5], [0.1,0.1], [0.1,0.1], [0.5,0.5] ] ])
Epsk = np.array([ [-1.5,-1.5], [-0.5,-0.5], [0.5,0.5], [1.5,1.5] ])


w0 = np.array([ 2.0, 1.5 ])
W  = np.array([ [0.6, 0.4], ])

ml = 1
Nk = 4
Ns = len(Eps0) + len(Epsk)

sectors = np.array([[0,0],])

data = h5py.File("input.h5", "w");

beta = data.create_dataset("BETA", shape=(), dtype='f', data=10.0)

hop_g = data.create_group("sectors")
hop_g.create_dataset("values", data=sectors)

bath = data.create_group("Bath")

if(Epsk.shape != Vk[0].shape):
    raise "Incorrect shape for Hybridisation and Epsk"
Epsk_g = bath.create_group("Epsk")
Epsk_g.create_dataset("values", shape=(len(Epsk),2,), data=Epsk, dtype=np.float)
Vk_g = bath.create_group("Vk")
Vk_g.create_dataset("values", data=np.array(Vk), dtype=np.float)

tk_g = data.create_group("tk")
tk_g.create_dataset("values", data=np.array(tk), dtype=np.float)

w0_g = bath.create_group("w0")
w0_g.create_dataset("values", shape=(len(w0),), data=w0, dtype=np.float)
W_g = bath.create_group("W")
W_g.create_dataset("values", data=np.array(W), dtype=np.float)


hop_g = data.create_group("Eps0")
hop_g.create_dataset("values", data=Eps0)

int_ds = data.create_dataset("U", shape=(), data=U)
int_ds = data.create_dataset("V", shape=(), data=V)

int_ds = data.create_dataset("AVG", shape=(), data=avg)

#int_g = data.create_group("chemical_potential")
int_ds = data.create_dataset("mu", shape=(), data=xmu)

