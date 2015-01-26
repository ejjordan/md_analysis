from mpl_toolkits.mplot3d import axes3d
import numpy as np
from scipy.interpolate import griddata
import matplotlib.pyplot as plt
import numpy.ma as ma
from numpy.random import uniform, seed
from matplotlib import cm
import re
import math
import os

class hiller(object):

    def __init__(self, protein_system, time, mins=None, maxs=None):
        self.protein_system=protein_system
        self.time=time
        self.mins=mins
        self.maxs=maxs

    def read_hills(self,stride,path=os.getcwd(),interpolate_points=100,
                   method='linear'):
        self.stride=stride
        self.interpolate_points=interpolate_points
        self.method=method
        self.path=path
        if self.stride==None:
            file="fes_{0}.dat".format(self.time)
        else:
            file="fes_{0}_{1}.dat".format(self.time,self.stride)
        FILE=open('{0}/{1}'.format(path,file))
        lines=FILE.readlines()
        FILE.close
        
        cv1_distance=[]
        cv2_distance=[]
        energy=[]
        for line in lines[10:]:
            if(re.match('.+',line)):
                vector = re.split('\s+',line)
                cv1_distance.append(float(vector[1]))
                cv2_distance.append(float(vector[2]))
                energy.append(float(vector[3]))

        self.cv1_distance=cv1_distance
        self.cv2_distance=cv2_distance
        self.energy=energy

        # define grid.
        self.cv1_grid = np.linspace(self.mins[0],self.maxs[0],self.interpolate_points)
        self.cv2_grid = np.linspace(self.mins[1],self.maxs[1],self.interpolate_points)
        # grid the data.
        self.energy_grid = griddata((self.cv1_distance, self.cv2_distance), self.energy,
                                    (self.cv1_grid[None,:], self.cv2_grid[:,None]),
                                    method=self.method)

    def draw_fig(self):
        # contour the gridded data, plotting dots at the randomly spaced data points.
        self.CS = plt.contour(self.cv1_grid,self.cv2_grid,self.energy_grid,10,
                              linewidths=0.5,colors='k')
        self.CS = plt.contourf(self.cv1_grid,self.cv2_grid,self.energy_grid,100,
                               cmap=plt.cm.jet)
        plt.colorbar() # draw colorbar
        plt.suptitle('{0} after {1} nanoseconds of metadynamics'
                     .format(self.protein_system, self.time))
        plt.xlabel('c-loop/a-loop distance (angstroms)')
        plt.ylabel('a-loop/aC helix distance (angstroms)')
        plt.show()

    def save_fig(self):
        # contour the gridded data, plotting dots at the randomly spaced data points.
        self.CS = plt.contour(self.cv1_grid,self.cv2_grid,self.energy_grid,10,
                              linewidths=0.5,colors='k')
        self.CS = plt.contourf(self.cv1_grid,self.cv2_grid,self.energy_grid,100,
                               cmap=plt.cm.jet)
        plt.colorbar() # draw colorbar
        if self.stride==None:
            plt.suptitle('{0} after {1} nanoseconds of metadynamics'
                         .format(self.protein_system, self.time))
        else:
            plt.suptitle('{0} after {1} nanoseconds of metadynamics'
                         .format(self.protein_system, (self.stride+1)))

        plt.xlabel('c-loop/a-loop distance (angstroms)')
        plt.ylabel('a-loop/aC helix distance (angstroms)')
        if self.stride==None:
            plt.savefig('{0}-{1}.png'.format(self.protein_system,self.time),format='png')
        else:
            plt.savefig('{0}-{1}_{2}.png'.format(self.protein_system,self.time,
                                                 self.stride),format='png')
        plt.clf()

    def read_hills_singleCV(self,stride=None, path=os.getcwd(),
                            interpolate_points=100, cv_index=1, energy_index=2,
                            single=True):
        self.stride=stride
        self.interpolate_points=interpolate_points
        self.path=path
        if self.stride==None:
            file="fes_{0}.dat".format(self.time)
        else:
            file="fes_{0}_{1}.dat".format(self.time,self.stride)
        FILE=open('{0}/{1}'.format(path,file))
        lines=FILE.readlines()
        FILE.close
        
        cv_distance=[]
        energy=[]
        for line in lines:
            if(re.match('#!',line)):
                pass
            else:
                if(re.match('.+',line)):
                    vector = re.split('\s+',line)
                    cv_distance.append(float(vector[cv_index]))
                    energy.append(float(vector[energy_index]))

        self.cv_distance=cv_distance
        self.energy=energy
        #plt.colorbar() # draw colorbar
        if single==True:
            line = plt.plot(self.cv_distance, self.energy)
            plt.suptitle('{0} after {1} nanoseconds of metadynamics'
                         .format(self.protein_system, self.time))
            plt.xlabel('CV')
            plt.ylabel('energy')
            plt.show()
        else:
            return self.cv_distance,self.energy

    def cv_progress(self, stride):
        data=[]
        for i in range(0,stride):
            cv,energy=self.read_hills_singleCV(stride=i,single=False)
            plt.plot(cv,energy, label="{0}".format(i))
            plt.legend()
        plt.suptitle('{0} after {1} nanoseconds of metadynamics'
                     .format(self.protein_system, self.time))
        plt.xlabel('CV')
        plt.ylabel('energy')
        plt.show()

    def frames(self, stride):
        #make one frame
        if stride == None:
            self.read_hills()
            self.save_fig()

        else:
            num_files=range(0,stride)
            for i in num_files:
                self.read_hills(stride=i)
                self.save_fig()
                print "wrote frame {0}".format(i)
                
