from MDAnalysis import *
import MDAnalysis
from MDAnalysis.analysis.align import *
import MDAnalysis.analysis.hbonds
import os
import sys
import Bio.PDB as bp
import numpy as np

"""
This script contains functions useful for analyzing MD trajectories
to determine if mutations (or other perturbations) cause a
structural change in a protein that may bias it towards a
different state than the wild type (or non-perturbed) system.

It is currently a work in progress!
"""


class analyzer(object):

    def __init__(self, topology, trajectory, path=os.getcwd()):
        self.path = path + '/'
        self.trajectory = trajectory
        self.topology = topology
        self.universe = MDAnalysis.Universe(self.path + self.topology,
                                            self.path + self.trajectory)

    def get_rms(self, ref_name="frame0.gro"):
        self.ref_name = ref_name
        gro_writer=MDAnalysis.coordinates.GRO.GROWriter(self.path
                                                        + self.ref_name)
        gro_writer.write(self.universe,0)
        self.ref_struct=Universe(self.path+self.ref_name)
        rms_fit_trj(self.universe,self.ref_struct,
                    {'mobile':'backbone', 'reference':'backbone'},
                    filename=None,rmsdfile=self.path + "fit.dat")
        os.remove(sys.path + "rmsfit_" + self.trajectory) #remove fitted traj


    def get_sasa(self, dssp_loc="/usr/bin/mkdssp"):
        self.dssp_loc = dssp_loc
        SASA={}
        protein=self.universe.selectAtoms("protein")
        for ts in self.universe.trajectory:
            writer=MDAnalysis.Writer("tmp.pdb")
            writer.write(protein)
            writer.close()
            parser=bp.PDBParser()
            structure=parser.get_structure('tmp','tmp.pdb')
            dssp=bp.DSSP(structure[0],'tmp.pdb',self.dssp_loc)
            for key in dssp.keys():
                resnum=dssp[key][0]
                sasa=dssp[key][2]
                if resnum.id[1] in SASA:
                    SASA[resnum.id[1]].append(sasa)
                else:
                    SASA[resnum.id[1]]=[sasa]
        count=0
        for key in SASA:
            print protein.resnames()[count],key,np.mean(SASA[key]),np.std(SASA[key],ddof=1)
            count+=1


    def get_salts(self):
        ions = self.universe.selectAtoms("type NA") + self.universe.selectAtoms("type CL")
        sel_basic = "(resname ARG or resname LYS) and (name NH* or name NZ)"
        sel_acidic = "(resname ASP or resname GLU) and (name OE* or name OD*)"
        basic=self.universe.selectAtoms(sel_basic)
        acidic=self.universe.selectAtoms(sel_acidic)
        #set up analysis of native contacts ("salt bridges"); with distance <6 A
        CA1 = MDAnalysis.analysis.contacts.ContactAnalysis1(self.universe, selection=(sel_acidic, sel_basic), refgroup=(acidic, basic), radius=3.0, outfile=self.path+"qsalt.dat")
        CA1.run(force=True)
        CA1.plot(self.path+"plot.pdf")


    def get_hbonds(self, distance=3.0, angle=160.0):
        sys.path.insert(0,"/home/joe/git/python_mdanalysis/mdanalysis/package/MDAnalysis/analysis")
        self.distance=distance
        self.angle=angle
        import joe_hbonds as joe
        hbonds=joe.HydrogenBondAnalysis(self.universe, 'protein', 'protein',
                                        update_selection1=True, 
                                        update_selection2=True, 
                                        detect_hydrogens='distance', 
                                        distance = self.distance,
                                        angle = self.angle, 
                                        distance_type="heavy")
        hbonds.run()
        hbonds_occupancy=hbonds.count_by_type()
        print hbonds_occupancy[0]
        hbonds_occupancy.generate_table()
        hbonds_occupancy.save_table(self.path+"hbonds-occupancy.pickle")

    def make_surface_graph(self, CaCa_dist=0.6, Ca1Cb2_dist=0.75,
                           Ca2Cb1_dist=0.75, dihedral=0.35):
        import scipy.spatial.distance as dist
        vector_list=[self.universe.selectAtoms("name CA").coordinates(),
                     self.universe.selectAtoms("name N and resname GLY or name CB").coordinates()]
        close_calphas=[]
        for calpha1_index in range(len(vector_list[0])):
            for calpha2_index in range(len(vector_list[0])):
                calpha1=vector_list[0][calpha1_index]
                calpha2=vector_list[0][calpha2_index]
                calphas_distance=dist.euclidean(calpha1,calpha2)
                if calphas_distance<10:
                    if calphas_distance==0:
                        continue
                    calpha1_resname=self.universe.selectAtoms("resid {0}".format(calpha1_index))
                    calpha2_resname=self.universe.selectAtoms("resid {0}".format(calpha1_index))
                    cbeta1=vector_list[1][calpha1_index]
                    cbeta2=vector_list[1][calpha2_index]
                    cbetas_distance=dist.euclidean(cbeta1,cbeta2)
                    calpha1_cbeta2_dist=dist.euclidean(calpha1,cbeta2)
                    calpha2_cbeta1_dist=dist.euclidean(calpha2,cbeta1)
                    #dihedral=MDAnalysis.core.distances.calc_torsions(calpha1,cbeta1,calpha2,cbeta2)
                    close_calphas.append([calpha1_index,calpha1_resname,
                                          calpha1,calpha2_index,
                                          calpha2_resname,  calpha2,cbeta1,
                                          cbeta2,calphas_distance,
                                          cbetas_distance,calpha1_cbeta2_dist,
                                          calpha2_cbeta1_dist])#,dihedral])
        self.calphas=vector_list[0]
        self.cbetas=vector_list[1]
        self.close_calphas=close_calphas
        
