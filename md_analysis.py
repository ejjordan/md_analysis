from MDAnalysis import *
import MDAnalysis
from MDAnalysis.analysis.align import *
import os

path="/home/joe/egfr/g719s/inactive-monomer/simulate/"
trj_file="md.part0001.xtc"
gro_file="md.part0001.gro"
frame0="frame0.gro"

trj=Universe(path+gro_file,path+trj_file)
gro_writer=MDAnalysis.coordinates.GRO.GROWriter(path+frame0)
gro_writer.write(trj,0)
ref=Universe(path+frame0)
rms_fit_trj(trj,ref,{'mobile':'backbone', 'reference':'backbone'},
            filename=None,rmsdfile=path+"fit.dat")

os.remove(path+"rmsfit_"+trj_file)

protein=trj.selectAtoms("protein")
ions = trj.selectAtoms("type NA") + trj.selectAtoms("type CL")
notwater=trj.selectAtoms("not resname SOL")



#MAY need these to take periodicity into account for atomselect
MDAnalysis.core.flags['use_periodic_selections'] = True
MDAnalysis.core.flags['use_KDTree_routines'] = False
