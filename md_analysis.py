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

    def __init__(self,topology,trajectory,start_frame=0,last_frame=None,
                 out_path=os.getcwd(),out_file=None):
        self.out_path = out_path + '/'
        self.trajectory = trajectory
        self.topology = topology
        self.start_frame=start_frame
        self.last_frame=last_frame
        self.universe = MDAnalysis.Universe(self.topology, self.trajectory)
        self.out_file=out_file
        self.universe.trajectory.start_timestep=self.start_frame
        if last_frame is not None:
            self.universe.trajectory.numframes=self.last_frame

    def get_rms(self, ref_name="frame0.gro"):
        self.ref_name = ref_name
        gro_writer=MDAnalysis.coordinates.GRO.GROWriter(self.out_path
                                                        + self.ref_name)
        gro_writer.write(self.universe,0)
        self.ref_struct=Universe(self.out_path+self.ref_name)
        rms_fit_trj(self.universe,self.ref_struct,
                    {'mobile':'backbone', 'reference':'backbone'},
                    filename=None, rmsdfile=self.out_path +
                    self.out_file + "_fit.dat")
        #os.remove(self.path + "rmsfit_" + self.trajectory) #remove fitted traj

    def get_sasa(self, dssp_loc="/usr/bin/mkdssp",skip=None):
        self.dssp_loc = dssp_loc
        SASA={}
        protein=self.universe.selectAtoms("protein")
        for ts in self.universe.trajectory:
            if skip:
                self.universe.trajectory.skip=skip
            sys.stdout.flush()
            sys.stdout.write('\rsasa [step {0}]  '.format(
                    self.universe.trajectory.frame))
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
        fp=open(self.out_path +self.out_file+ "_sasa.dat",'w')
        for key in SASA:
            fp.write("{0}\t{1}\t{2}\t{3}\n".format(protein.resnames()[count],key,np.mean(SASA[key]),np.std(SASA[key],ddof=1)))
            count+=1
        fp.close()
        sys.stdout.write('\rSASA table created     ')
        print


    def get_salts_lame(self):
        """
        This version of the get salts routine uses the MDAnalysis contact
        analysis routine. It is functional but the output is clunky so it
        should only be used if you want to learn how to parse this output.
        """
        ions = self.universe.selectAtoms("type NA") + self.universe.selectAtoms("type CL")
        sel_basic = "(resname ARG or resname LYS) and (name NH* or name NZ)"
        sel_acidic = "(resname ASP or resname GLU) and (name OE* or name OD*)"
        basic=self.universe.selectAtoms(sel_basic)
        acidic=self.universe.selectAtoms(sel_acidic)
        #set up analysis of native contacts ("salt bridges"); with distance <6 A
        CA1 = MDAnalysis.analysis.contacts.ContactAnalysis1(self.universe, selection=(sel_acidic, sel_basic), refgroup=(acidic, basic), radius=3.0, outfile=self.out_path+"qsalt.dat")
        CA1.run(force=True)
        return CA1
        #CA1.plot(self.out_path+"plot.pdf")

    def get_salts(self, distance=4.0, angle=60.0):
        """
        This version of the get salts routine usese a modified version of the
        MDAnalysis hbonds routine. The forcefield type used is 'ion', which
        includes nitrogens and oxygens in the side chain of HIS, GLN, GLU, 
        ASP, ARG, and LYS. Note that currently salt bridges to solvent ions
        will not be found.
        """
        sys.path.insert(0,"/home/joe/git/md_analysis/")
        self.salt_distance=distance
        self.salt_angle=angle
        import joe_hbonds as joe
        salts=joe.HydrogenBondAnalysis(self.universe, 'protein', 'protein',
                                       update_selection1=True, 
                                       update_selection2=True, 
                                       detect_hydrogens='distance', 
                                       distance = self.salt_distance,
                                       angle = self.salt_angle, 
                                       distance_type="heavy",
                                       forcefield='ion')
        salts.run()
        salts_occupancy=salts.count_by_type()
        salts_occupancy_list=salts_occupancy.tolist()
        fp=open(self.out_path +self.out_file+ "_salts.dat",'w')
        for i in range(len(salts_occupancy_list)):
            fp.write("{0}\n".format(salts_occupancy_list[i]))
        fp.close()

    def get_hbonds(self, distance=3.0, angle=160.0):
        sys.path.insert(0,"/home/joe/git/md_analysis/")
        import joe_hbonds as joe
        self.hbond_distance=distance
        self.hbond_angle=angle
        hbonds=joe.HydrogenBondAnalysis(self.universe, 'protein', 'protein',
                                        update_selection1=True, 
                                        update_selection2=True, 
                                        detect_hydrogens='distance', 
                                        distance = self.hbond_distance,
                                        angle = self.hbond_angle, 
                                        distance_type="heavy")
        hbonds.run()
        hbonds_occupancy=hbonds.count_by_type()
        hbonds_occupancy_list=hbonds_occupancy.tolist()
        fp=open(self.out_path +self.out_file+ "_hbonds.dat",'w')
        for i in range(len(hbonds_occupancy_list)):
            fp.write("{0}\n".format(hbonds_occupancy_list[i]))
        fp.close()
