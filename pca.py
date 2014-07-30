from scipy import linalg as la
from MDAnalysis import *
import MDAnalysis
import MDAnalysis.analysis as analysis
import numpy as np
import os
import matplotlib.pylab as plt

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
        ref="ref.gro"
        gro_writer=MDAnalysis.coordinates.GRO.GROWriter(self.targetdir+ref)
        gro_writer.write(self.universe,0)
        reference=Universe(self.targetdir+ref).selectAtoms(self.selection)
        num_atoms=selection.numberOfAtoms()
        num_frames=self.universe.trajectory.numframes
        dof=num_atoms*3
        cov=np.zeros((dof,dof))
        
        print "covariancs matrix will have shape {0} from {1} frames".format(np.shape(cov),num_frames)
        
        coordsum=np.zeros(dof)
        avg_struc=np.zeros(num_atoms,3)
        for ts in self.universe.trajectory:
            analysis.align.alignto(selection,reference,mass_weighted=True)
            coords = selection.coordinates().flatten()
            coordsum += coords
            cov += np.outer(coords,coords)
            avg_struc += selection.coordinates()
        cov /= num_frames            
        coordsum /= num_frames
        avg_struc /= num_frames
        self.average_structure=avg_struc
        cov -= np.outer(coordsum,coordsum)
        self.covariance=cov
        self.eigvals,self.eigvecs=la.eig(self.covariance)

    def get_covariance(self):
        return self.covariance

    def get_eigenvalues(self):
        return self.eigvals

    def get_eigenvectors(self):
        return self.eigvecs

    def get_average_structure(self):
        return self.average_structure

"""
    def make_PCA_traj(self, pca_file="PCA.xtc"):
        from xdrfile.TRR import TRRReader, TRRWriter, Timestep

        # 'modes' is a mode object with M PCs, similar to a MxNx3 array
        # 'xav' the average coordinates, a Nx3 array for N atoms

        N = len(self.average_structure)   # number of atoms

        W = Writer(self.targetdir + pca_file, numatoms=N)            # TRR writer
        ts = MDAnalysis.coordinates.TRR.Timestep(N)  # TRR time step

        for frame,mode in enumerate(self.eigvecs[0:9]):
            ts.lmbda = -1
            if frame<=1:
                ts._pos[:] = xav
            else:
                ts._pos[:] = mode.scaledToNorm(1.).array*10   # nm to angstroms
                ts.frame = frame         # manually change the frame number
                ts.step = frame - 1
                if frame <= 1:
                    ts.time = frame-1
                else:
                    ts.time = mode.frequency
                    W.write(ts)             # converts angstrom to nm for gmx

        W.close()
"""
