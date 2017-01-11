#!/usr/bin/env python

import h5py
import numpy as np

import U_matrix

U = np.real(U_matrix.U_matrix(2, U_int=6.6, J_hund=0.9, basis='cubic'))
xmu = 44.44
Eps0 = np.array([[-0.701786,-0.701786], [-0.812433,-0.812433], [-0.573927,-0.573927], [-0.812433,-0.812433], [-0.701786,-0.701786] ])
Vk =   [ np.array([ [0.40219,  0.40219], [ 0.369274,0.369274]]),
         np.array([ [0.731886,0.731886], [ 0.480395,0.480395]]),
         np.array([ [0.484533,0.484533], [ 0.293002,0.293002]]),
         np.array([ [0.731886,0.731886], [ 0.480395,0.480395]]),
         np.array([ [0.40219,  0.40219], [ 0.369274,0.369274]]) ]
Epsk = [ np.array([ [-1.74659,-1.74659], [-0.937719,-0.937719]]),
         np.array([ [-1.73124,-1.73124], [ -1.11657,-1.11657] ]),
         np.array([ [-1.93049,-1.93049], [-1.20792,-1.20792]]),
         np.array([ [-1.73124,-1.73124], [-1.11657,-1.11657]]),
         np.array([ [-1.74659,-1.74659], [-0.937719,-0.937719]] ) ]

Ns = 15
ml = 5
sectors = np.array([[14,14],])

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

hop_g = data.create_group("Eps0")
hop_g.create_dataset("values", data=Eps0)

int_g = data.create_group("interaction")
int_ds = int_g.create_dataset("values", shape=(ml,ml,ml,ml,), data=U)

#int_g = data.create_group("chemical_potential")
int_ds = data.create_dataset("mu", shape=(), data=xmu)

