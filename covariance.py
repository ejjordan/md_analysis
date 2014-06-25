
from scipy import linalg as la
from MDAnalysis import *
import numpy as np

path="/home/joe/egfr/g719s/inactive-monomer/simulate/"
trj_file="short.xtc"
gro_file="md.part0001.gro"
system=Universe(path+gro_file,path+trj_file)

#system = AtomGroup.System("data/beta_op_adp-vmd.psf", "data/step1.dcd")
asel = system.selectAtoms('protein and backbone and not type O')

masses = np.repeat(asel.masses(), 3)
mass_matrix = np.sqrt(np.identity(len(masses))*masses)

skip = 1
num_ts = system.trajectory.numframes/skip
num_coor = len(asel)*3

ca_pos = system.trajectory.timeseries(asel, skip=skip, format='fac')
collection.addTimeseries(Timeseries.CenterOfMass(asel))
ca_com = np.transpose(system.trajectory.correl(collection, skip=skip))

# Recenter on center of mass of selection
ca_pos -= ca_com[:,np.NewAxis]

# Remove rotational degrees of freedom
ref = ca_pos[0]
for coor in ca_pos[1:]:
    rotmatrix = rms_fitting.rms_rotation_matrix(coor, ref, kg)
    coor[:] = np.matrixmultiply(coor, rotmatrix).astype(np.Float32)

ca = np.reshape(ca_pos, (num_ts, -1))
ca_avg = np.average(ca)
ca2 = ca - ca_avg[np.NewAxis,:]
ca_cov = np.zeros((num_coor, num_coor), np.Float)
for ts in ca2:
    ca_cov += np.outerproduct(ts, ts)
ca_cov /= num_ts
ca_cov1 = np.matrixmultiply(ca_cov, mass_matrix)
del ca_cov
ca_cov2 = np.matrixmultiply(mass_matrix, ca_cov1)
del ca_cov1

N_av = 6.0221367e23
hplanck_bar = 6.6260755e-34/(2*np.pi)
k =  1.3806580000000001e-23
T = 300 # kelvin
eigenv, eigenvec = la.eigenvectors(ca_cov2)
real = [e.real/100. for e in eigenv]
f = open('eigenval.dat', 'w')
for i, val in enumerate(real):
    f.write(`i+1` + '\t' + `val` + '\n')
f.close()

eigenval = eigenv*1.6605402e-27*1e-20
omega_i = np.sqrt(k*T/(eigenval))

term = (hplanck_bar*omega_i)/(k*T)
summation_terms = (term/(np.exp(term)-1.))-np.log(1.-np.exp(-term))
S_ho = k*N_av*np.sum(summation_terms)
