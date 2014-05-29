from MDAnalysis import *
u = Universe("md.part0001.gro", "md.part0001.xtc")

#MAY need these to take periodicity into account for atomselect
MDAnalysis.core.flags['use_periodic_selections'] = True
MDAnalysis.core.flags['use_KDTree_routines'] = False
