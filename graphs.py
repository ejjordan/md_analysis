import MDAnalysis
import numpy as np
from pygraph.classes.graph import graph
from pygraph.readwrite.dot import write
import MDAnalysis.KDTree.NeighborSearch as ns

class graphs(object):

    def __init__(self, universe):
        self.universe=universe

    def calphas(self,distance=6):
        self.ca_graph=graph()
        self.distance=distance
        ca=self.universe.selectAtoms("name CA")
        ca_coords=ca.coordinates()
        ca_neighbors=ns.CoordinateNeighborSearch(ca_coords)
        for resid in ca.resids():
            self.ca_graph.add_node("{0}".format(resid))
            self.ca_graph.add_node_attribute("{0}".format(resid),
                                        ca.resnames()[resid])

        for resid in ca.resids():
            neighbors=ca_neighbors.search(ca_coords[resid],self.distance,
                                          distances=True)
            for index in range(len(neighbors[0])):
                if resid==neighbors[0][index] or resid==neighbors[0][index]+1 or resid==neighbors[0][index]-1:
                    continue
                try:
                    self.ca_graph.add_edge(("{0}".format(resid),"{0}".format(neighbors[0][index])),wt=neighbors[1][index])
                except:
                    pass

    def draw_graph(self,filename="protein.png"):
        import gv
        dot = write(self.ca_graph)
        gvv = gv.readstring(dot)
        gv.layout(gvv,'dot')
        gv.render(gvv,'png',filename)
