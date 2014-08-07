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

    def dagget_residue(self,start_res=0,C_distance=5.4,non_C_distance=4.6):
        self.residue_graph=igraph.Graph()
        self.C_distance=C_distance
        self.non_C_distance=non_C_distance
        protein=self.universe.selectAtoms("protein")
        carbons=self.universe.selectAtoms("protein and name C*")
        for resid in protein.resids():
            residue=self.universe.selectAtoms("resid {0}".format(resid))
            self.residue_graph.add_vertex("{0} {1}".format(
                    residue.resnames()[0],resid+start_res))

        search=ns.CoordinateNeighborSearch(protein.coordinates())
        neighbors=search.search_all(self.non_C_distance)
        for index in range(len(neighbors)):
            atom1,atom2=neighbors[index]
            atom1_type=protein[atom1].name
            atom1_resnum=protein[atom1].resid
            atom2_type=protein[atom2].name
            atom2_resnum=protein[atom2].resid
#            if atom1_resnum == "C" and atom2_resnum == "C"
#                continue
            if atom1_resnum == atom2_resnum:
                continue
            self.residue_graph.add_edge(atom1_resnum,atom2_resnum)

        search=ns.CoordinateNeighborSearch(carbons.coordinates())
        neighbors=search.search_all(self.C_distance)
        for index in range(len(neighbors)):
            atom1,atom2=neighbors[index]
            atom1_type=protein[atom1].name
            atom1_resnum=protein[atom1].resid
            atom2_type=protein[atom2].name
            atom2_resnum=protein[atom2].resid
            if atom1_resnum == atom2_resnum:
                continue
            self.residue_graph.add_edge(atom1_resnum,atom2_resnum)

    def draw_cricle(self):
        color_dict = {"ach": "red", "aloop": "green", "ploop": "orange",
                      "cloop": "pink", "None": "blue"}
        self.residue_graph.vs[55:75]["region"]='ach'
        ach_edges=self.residue_graph.es.select(_source_in=
                                               self.residue_graph.vs[55:75])
        ach_edges["color"]="red"
        ach_edges2=self.residue_graph.es.select(_source_in=
                                               self.residue_graph.vs[55:75])
        ach_edges["color"]="red"
        self.residue_graph.vs[155:185]["region"]='aloop'
        aloop_edges=self.residue_graph.es.select(_source_in=
                                                 self.residue_graph.vs[155:185])
        aloop_edges["color"]="green"
        self.residue_graph.vs[23:30]["region"]='ploop'
        ploop_edges=self.residue_graph.es.select(_source_in=
                                                 self.residue_graph.vs[23:30])
        ploop_edges["color"]="orange"
        self.residue_graph.vs[136:147]["region"]='cloop'
        cloop_edges=self.residue_graph.es.select(_source_in=
                                               self.residue_graph.vs[136:147])
        cloop_edges["color"]="pink"
        for index in range(self.residue_graph.vcount()):
            self.residue_graph.vs[index]["color"]=color_dict["{0}".format(
                    self.residue_graph.vs[index]["region"])]
        layout=self.residue_graph.layout_circle()
        igraph.plot(self.residue_graph,layout=layout)
