from MDAnalysis import *
import MDAnalysis
from MDAnalysis.analysis.align import *
import MDAnalysis.analysis.hbonds
import os
import sys

#MAY need these to take periodicity into account for atomselect
#MDAnalysis.core.flags['use_periodic_selections'] = True
#MDAnalysis.core.flags['use_KDTree_routines'] = False

path="/home/joe/egfr/g719s/inactive-monomer/simulate/"
trj_file="short.xtc"
gro_file="md.part0001.gro"
frame0="frame0.gro"

trj=Universe(path+gro_file,path+trj_file)
gro_writer=MDAnalysis.coordinates.GRO.GROWriter(path+frame0)
gro_writer.write(trj,0)
ref=Universe(path+frame0)
#rms_fit_trj(trj,ref,{'mobile':'backbone', 'reference':'backbone'},
#            filename=None,rmsdfile=path+"fit.dat")
#os.remove(path+"rmsfit_"+trj_file)

protein=trj.selectAtoms("protein")
ions = trj.selectAtoms("type NA") + trj.selectAtoms("type CL")
notwater=trj.selectAtoms("not resname SOL")

sys.path.insert(0,"/home/joe/git/python_mdanalysis/mdanalysis/package/MDAnalysis/analysis")

SASA={}
import Bio.PDB as bp
import numpy as np
for ts in trj.trajectory:
    writer=MDAnalysis.Writer("tmp.pdb")
    writer.write(protein)
    writer.close()
    parser=bp.PDBParser()
    structure=parser.get_structure('tmp','tmp.pdb')
    dssp=bp.DSSP(structure[0],'tmp.pdb',"/usr/bin/mkdssp")
    for key in dssp.keys():
#        print key
        resnum=dssp[key][0]
        sasa=dssp[key][2]
        #print resnum.id[1],sasa
        if resnum.id[1] in SASA:
            SASA[resnum.id[1]].append(sasa)
        else:
            SASA[resnum.id[1]]=[sasa]
#print np.mean(SASA[697])
#print np.std(SASA[697],ddof=1)
count=0
for key in SASA:
    print protein.resnames()[count],key,np.mean(SASA[key]),np.std(SASA[key],ddof=1)
    count+=1

"""
import joe_hbonds as joe
hbonds=joe.HydrogenBondAnalysis(trj, 'protein', 'protein',
                                update_selection1=True, 
                                update_selection2=True, 
                                detect_hydrogens='distance', 
                                distance=3.00, angle=160.0, 
                                distance_type="heavy")
hbonds.run()
hbonds_occupancy=hbonds.count_by_type()
print hbonds_occupancy[0]
#hbonds_occupancy.generate_table()
#hbonds_occupancy.save_table(path+"hbonds-occupancy.pickle")


sel_basic = "(resname ARG or resname LYS) and (name NH* or name NZ)"
sel_acidic = "(resname ASP or resname GLU) and (name OE* or name OD*)"
basic=trj.selectAtoms(sel_basic)
acidic=trj.selectAtoms(sel_acidic)
# set up analysis of native contacts ("salt bridges"); with distance <6 A
CA1 = MDAnalysis.analysis.contacts.ContactAnalysis1(trj, selection=(sel_acidic, sel_basic), refgroup=(acidic, basic), radius=3.0, outfile=path+"qsalt.dat")
CA1.run(force=True)
CA1.plot(path+"plot.pdf")
"""
