from scipy import linalg as la
from MDAnalysis import *
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
        trajectorybase = os.path.splitext(os.path.basename(trajectory))[0]
        output = trajectorybase + '_pca.dat'
        self.output = os.path.join(self.targetdir, output)
        self.universe = Universe(topology, trajectory)
    
    def make_covariance(self):

        selection = self.universe.selectAtoms(self.selection)
        num_atoms=selection.numberOfAtoms()
        num_frames=self.universe.trajectory.numframes
        dof=num_atoms*3
        cov=np.zeros((dof,dof))
        
        print "covariancs matrix will have shape {0} from {1} frames".format(np.shape(cov),num_frames)
        
        num_confs=0
        coordsum=np.zeros(dof)
        for ts in self.universe.trajectory:
            frame = ts.frame
            coords=selection.coordinates().flatten()
            coordsum += coords
            cov +=np.outer(coords,coords)
            num_confs+=1
        cov /= num_confs
        coordsum /= num_confs
        mean = coordsum
        cov -= np.outer(coordsum,coordsum)
        print cov

