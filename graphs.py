import MDAnalysis
import numpy as np
from pygraph.classes.graph import graph
from pygraph.readwrite.dot import write
import MDAnalysis.KDTree.NeighborSearch as ns
import igraph

class graphs(object):

    def __init__(self, universe):
        self.universe=universe


    def calphas(self,distance=6,start_res=0,take_connected=False):
        self.ca_graph=igraph.Graph()
        self.distance=distance
        self.take_connected=take_connected
        ca=self.universe.selectAtoms("name CA")
        ca_coords=ca.coordinates()
        ca_neighbors=ns.CoordinateNeighborSearch(ca_coords)
        for resid in ca.resids():
            self.ca_graph.add_vertex("{0} {1}".format(ca.resnames()[resid],
                                                      resid+start_res))

        if self.take_connected==False:
            for resid in ca.resids():
                neighbors=ca_neighbors.search(ca_coords[resid],self.distance,
                                              distances=True)
                for index in range(len(neighbors[0])):
                    if      resid==neighbors[0][index] or \
                            resid==neighbors[0][index]+1 or \
                            resid==neighbors[0][index]-1:
                        continue
                    self.ca_graph.add_edge(resid,neighbors[0][index],
                                           weight=neighbors[1][index])

        if self.take_connected==True:
            for resid in ca.resids():
                neighbors=ca_neighbors.search(ca_coords[resid],self.distance,
                                              distances=True)
                for index in range(len(neighbors[0])):
                    if resid==neighbors[0][index]:
                        continue
                    self.ca_graph.add_edge(resid,neighbors[0][index],
                                           weight=neighbors[1][index])

    def dagget_residue(self,start_res=0):
        self.residue_graph=igraph.Graph()
        protein=self.universe.selectAtoms("protein")
        for resid in protein.resids():
            residue=self.universe.selectAtoms("resid {0}".format(resid))
            self.residue_graph.add_vertex("{0} {1}".format(
                    residue.resnames()[0],resid+start_res))

        for resid1 in protein.resids():
            residue1=self.universe.selectAtoms("resid {0}".format(resid1))
            for resid2 in protein.resids():
                residue2=self.universe.selectAtoms("resid {0}".format(resid2))
                

            res_atoms=residue.atoms()
            res_neighbors=ns.CoordinateNeighborSearch(res_coords)
            for resid in ca.resids():
            neighbors=ca_neighbors.search(ca_coords[resid],self.distance,
                                          distances=True)
            for index in range(len(neighbors[0])):
                if resid==neighbors[0][index]:
                    continue
                self.ca_graph.add_edge(resid,neighbors[0][index],
                                       weight=neighbors[1][index])
