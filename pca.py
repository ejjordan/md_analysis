from scipy import linalg as la
from MDAnalysis import *
import MDAnalysis
import MDAnalysis.analysis as analysis
import numpy as np
import os

class PCA(object):
    """perform principal component analysis"""
    def __init__(self, topology, trajectory, selection='protein and backbone', 
                 targetdir=os.path.curdir):
        self.topology = topology
        self.trajectory = trajectory
        self.targetdir = targetdir
        self.selection = selection
        trajbase = os.path.splitext(os.path.basename(trajectory))[0]
        output = trajbase + '_pca.dat'
        self.output = os.path.join(self.targetdir, output)
        self.universe = Universe(topology, trajectory)
    
    def make_covariance(self):
        selection = self.universe.selectAtoms(self.selection)
        #reference = self.universe.selectAtoms(self.selection)
        ref="ref.gro"
        gro_writer=MDAnalysis.coordinates.GRO.GROWriter(self.targetdir+ref)
        gro_writer.write(self.universe,0)
        reference=Universe(self.targetdir+ref).selectAtoms(self.selection)
        num_atoms=selection.numberOfAtoms()
        num_frames=self.universe.trajectory.numframes
        dof=num_atoms*3
        cov=np.zeros((dof,dof))
        
        print "covariancs matrix will have shape {0} from {1} frames".format(np.shape(cov),num_frames)
        
        num_confs=0
        coordsum=np.zeros(dof)
        for ts in self.universe.trajectory:
            R=analysis.align.rotation_matrix(selection.coordinates(),reference.coordinates())
            selection.atoms.rotate(R[0])
            coords=selection.coordinates().flatten()
            coordsum += coords
            cov +=np.outer(coords,coords)
            num_confs+=1
        cov /= num_confs
        coordsum /= num_confs
        mean = coordsum
        cov -= np.outer(coordsum,coordsum)
        print cov

