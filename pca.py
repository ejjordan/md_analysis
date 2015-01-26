from scipy import linalg as la
from MDAnalysis import *
import MDAnalysis
import MDAnalysis.analysis as analysis
import numpy as np
import os,sys
import matplotlib.pylab as plt
import gzip

class PCA(object):
    """perform principal component analysis"""
    def __init__(self, topology, trajectory, align_selection='name CA',
                 pca_selection="name CA", start_frame=0, last_frame=None, 
                 targetdir=os.path.curdir, skip=None, outname=None):
        self.topology = topology
        self.trajectory = trajectory
        self.targetdir = targetdir
        self.skip = skip
        #trajbase = os.path.splitext(os.path.basename(trajectory))[0]
        #output = trajbase + '_pca.dat'
        self.output = self.targetdir + "/" + outname
        self.start_frame=start_frame
        self.last_frame=last_frame
        self.universe = Universe(topology, trajectory)
        self.universe.trajectory.start_timestep=self.start_frame
        if last_frame is not None:
            self.universe.trajectory.numframes=self.last_frame
        self.align_selection = align_selection
        self.pca_selection = self.universe.selectAtoms(pca_selection)

    
    def make_covariance(self):
        align_selection = self.universe.selectAtoms(self.align_selection)
        ref="ref.gro"
        gro_writer=MDAnalysis.coordinates.GRO.GROWriter(self.targetdir+ref)
        gro_writer.write(self.universe,0)
        reference=Universe(self.targetdir+ref).selectAtoms(self.align_selection)
        num_atoms=self.pca_selection.numberOfAtoms()
        num_frames=self.universe.trajectory.numframes
        dof=num_atoms*3
        cov=np.zeros((dof,dof))
        print "covariancs matrix will have shape {0} from {1} frames".format(np.shape(cov),num_frames)
        coordsum=np.zeros(dof)
        avg_struc=np.zeros((num_atoms,3))
        for ts in self.universe.trajectory:
            frame_num=self.universe.trajectory.frame
            if (frame_num % 10) == 0:
                sys.stdout.flush()
                sys.stdout.write('\rPCA [step {0}]  '.format(frame_num))
            if self.skip:
                self.universe.trajectory.skip=self.skip
            analysis.align.alignto(align_selection,reference,mass_weighted=True)
            coords = self.pca_selection.coordinates().flatten()
            coordsum += coords
            cov += np.outer(coords,coords)
            avg_struc += self.pca_selection.coordinates()
        cov /= num_frames            
        coordsum /= num_frames
        avg_struc /= num_frames
        self.average_structure=avg_struc
        cov -= np.outer(coordsum,coordsum)
        self.covariance=cov
        print "\ncomputing eigenvectors and eigenvalues"
        self.eigvals,self.eigvecs=la.eig(np.array(self.covariance))

    def write_data(self):
        fp=open("{0}_eigenvals.dat".format(self.output), 'w')
        for i in range(len(self.eigvals)):
            fp.write("{0}\n".format(self.eigvals[i]))
        fp.close()
        fp=gzip.GzipFile("{0}_eigenvecs.gz".format(self.output), 'wb')
        for i in range(len(self.eigvecs)):
            for j in range(len(self.eigvecs[i])):
                fp.write("{0} ".format(self.eigvecs[i][j]))
            fp.write("\n")
        fp.close()



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
