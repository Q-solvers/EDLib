import h5py
import numpy as np
import argparse


parser = argparse.ArgumentParser(description='Generate Anderson model input for ED')
parser.add_argument("--ns", type=int, help="Number of spin states", default=2)
parser.add_argument("--norb", type=int, help="Number of orbitals", default=1)
parser.add_argument("--nbath", type=int, help="Number of bath states", default=1)
parser.add_argument("--t", type=str, help="t.txt file")
parser.add_argument("--v", type=str, help="v.txt file")
parser.add_argument("--bath", type=str, help="bath.txt file")
parser.add_argument("--input_name", type=str, help="input.h5 file", default="input.h5")
parser.add_argument("--mu", type=float, help="chemical potential", default=0.0)
parser.add_argument("--factor", type=float, help="Interaction prefactor", default=1.0)
args = parser.parse_args()
assert(args.ns == 2 or args.ns == 1) # only 1 or 2 spin states are allowed


# read interaction
U = np.zeros((args.norb,args.norb,args.norb,args.norb), dtype=float)
with open(args.v) as f:
    lns = f.readlines()
    for ln in lns:      
        lnspl = ln.split()
        i = int(lnspl[0])
        j = int(lnspl[1])
        k = int(lnspl[3])
        l = int(lnspl[2])
        U[i,j,k,l] = float(lnspl[4]) * args.factor

# read hopping
H0 = np.zeros((args.norb,args.norb,args.ns))
with open(args.t) as f:
    lns = f.readlines()
    for ln in lns:
        lnspl = ln.split()
        i = int(lnspl[0])
        j = int(lnspl[1])
        t = float(lnspl[2])
        H0[i,j,0] = t
        H0[i,j,1] = t


# read bath
with open(args.bath) as f:
    lns = f.readlines()
    assert len(lns[0].split()) == args.nbath
    assert(len(lns) == args.norb + 1)

Epsk = np.zeros((args.nbath, args.ns))
Vk = np.zeros((args.norb, args.nbath, args.ns))

for ib, b in enumerate(lns[0].split()):
    Epsk[ib, :] = float(b)

for im, ln in enumerate(lns[1:]):
    lnspl = ln.split()
    for ib, b in enumerate(lnspl):
        Vk[im, ib, :] = float(b)


# write to hdf5
input_file = args.input_name
data = h5py.File(input_file, "w")

## write bath
bath = data.create_group("Bath")
bath["Epsk/values"] = Epsk

for i in range(args.norb):
    if(Epsk.shape != Vk[i].shape):
        raise "Incorrect shape for Hybridisation and Epsk"
    ## Write bath hybridization
    Vk_g = bath.create_group("Vk_" + str(i))
    Vk_g["values"] = Vk[i]
    ## write impurity noninteracting Hamiltonian
    t0_g = data.create_group("H0_" + str(i))
    t0_g["values"] = H0[i]

## write interaction
UU = np.zeros((args.ns,args.ns) + U.shape)
UU[:,:] = U

int_g = data.create_group("interaction")
int_ds = int_g.create_dataset("values", shape=(args.ns,args.ns,args.norb,args.norb,args.norb,args.norb,), data=UU)

## write chemical potential
mu_ds = data.create_dataset("mu", shape=(), data=args.mu)

## write orbital pairs
orb_pairs = np.array([[i,j] for i in range(args.norb) for j in range(args.norb)])
data["GreensFunction_orbitals/values"] = orb_pairs

