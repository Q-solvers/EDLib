#!/usr/bin/env python

import h5py
import numpy as np

U = np.array([6.6])
xmu = np.array([44.44])
Vk = np.array([0.56434, 0.68392, 0.29519])
epsk = np.array([-2.37325, -0.87328, 2.01265])

sectors = np.array([[2,2]])


Ns = len(U)
Nk = len(Vk)

data = h5py.File("input.h5", "w");

hop_g = data.create_group("sectors")
hop_g.create_dataset("values", data=sectors)


int_g = data.create_group("interaction")
int_ds = int_g.create_dataset("values", shape=(Ns,), data=U)

int_g = data.create_group("chemical_potential")
int_ds = int_g.create_dataset("values", shape=(Ns,), data=xmu)

int_g = data.create_group("Bath/Vk_1")
int_ds = int_g.create_dataset("values", shape=(Nk,), data=Vk)

int_g = data.create_group("Bath/Epsk_1")
int_ds = int_g.create_dataset("values", shape=(Nk,), data=epsk)
