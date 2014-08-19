import MDAnalysis
import numpy as np
from pygraph.classes.graph import graph
from pygraph.readwrite.dot import write
import MDAnalysis.KDTree.NeighborSearch as ns
import igraph

class graphs(object):

    def __init__(self, universe):
        """
        To initialize a graph object, provide an MDAnalysis universe object.
        You have several choices where to go from here. Different types of
        graph can be constructed by issuing the appropriate function call.
        Calling a graph of any type will result in a graph instance called
        'graff' for all graph types, for ease of function reuse.
        """
        self.universe=universe
        

    def calphas(self,distance=6,start_res=0,take_connected=False):
        self.graff=igraph.Graph()
        self.distance=distance
        self.take_connected=take_connected
        ca=self.universe.selectAtoms("name CA")
        ca_coords=ca.coordinates()
        ca_neighbors=ns.CoordinateNeighborSearch(ca_coords)
        for resid in ca.resids():
            self.graff.add_vertex("{0} {1}".format(ca.resnames()[resid],
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
                    self.graff.add_edge(resid,neighbors[0][index],
                                           weight=neighbors[1][index])

        if self.take_connected==True:
            for resid in ca.resids():
                neighbors=ca_neighbors.search(ca_coords[resid],self.distance,
                                              distances=True)
                for index in range(len(neighbors[0])):
                    if resid==neighbors[0][index]:
                        continue
                    self.graff.add_edge(resid,neighbors[0][index],
                                           weight=neighbors[1][index])

    def dagget_residue(self,start_res=0,C_distance=5.4,non_C_distance=4.6):
        self.graff=igraph.Graph()
        self.C_distance=C_distance
        self.non_C_distance=non_C_distance
        protein=self.universe.selectAtoms("protein")
        carbons=self.universe.selectAtoms("protein and name C*")
        for resid in protein.resids():
            residue=self.universe.selectAtoms("resid {0}".format(resid))
            self.graff.add_vertex("{0} {1}".format(
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
            eid=self.graff.get_eid(atom1_resnum,atom2_resnum,
                                           error=False)
            if eid < 0:
                self.graff.add_edge(atom1_resnum,atom2_resnum,weight=1)
            else:
                self.graff.es[eid]["weight"]+=1


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
            eid=self.graff.get_eid(atom1_resnum,atom2_resnum,
                                           error=False)
            if eid < 0:
                self.graff.add_edge(atom1_resnum,atom2_resnum,weight=1)
            else:
                self.graff.es[eid]["weight"]+=1


    def compute_graph_properties(self,graff):
        ceb=graff.community_edge_betweenness()
        clustering=ceb.as_clustering()
        modularity=clustering.modularity
        graff["modularity"]=modularity
        graff["avg_path_len"]=graff.average_path_length()
        degree_distribution=graff.degree_distribution()
        graff["degree distribribution mean"]=degree_distribution.mean
        graff["degree distribribution sd"]=degree_distribution.sd
        graff["density"]=graff.density(loops=True)
        graff["diamter"]=graff.diameter()
        graff["edge count"]=graff.ecount()
        laplacian=graff.laplacian()
        eigenvals=np.linalg.eigvals(laplacian)
        graff["laplacian"]=graff.laplacian()
        graff.vs["eigenvals"]=eigenvals
        graff["maxdegree"]=graff.maxdegree()
        graff["transitivity"]=graff.transitivity_avglocal_undirected(
            weights=graff.es["weight"])
        

    def draw_cricle(self):
        color_dict = {"ach": "red", "aloop": "green", "ploop": "blue",
                      "cloop": "orange", "None": "white"}
        all_verts=self.graff.vs
        ach_verts=self.graff.vs[55:75]
        aloop_verts=self.graff.vs[157:179]
        ploop_verts=self.graff.vs[23:30]
        cloop_verts=self.graff.vs[136:147]
        ach_edges=self.graff.es.select(_between=(ach_verts,all_verts))
        ach_edges["color"]="red"
        ach_verts["region"]='ach'
        aloop_edges=self.graff.es.select(_between=(aloop_verts,all_verts))
        aloop_edges["color"]="green"
        aloop_verts["region"]='aloop'
        ploop_verts["region"]='ploop'
        ploop_edges=self.graff.es.select(_between=(ploop_verts,all_verts))
        ploop_edges["color"]="blue"
        cloop_verts["region"]='cloop'
        cloop_edges=self.graff.es.select(_between=(cloop_verts,all_verts))
        cloop_edges["color"]="orange"
        for index in range(self.graff.vcount()):
            self.graff.vs[index]["color"]=color_dict["{0}".format(
                    self.graff.vs[index]["region"])]
        layout=self.graff.layout_circle()
        igraph.plot(self.graff,layout=layout)

    def draw_single_node(self, node):
        subgraph=self.graff.subgraph(self.neighborhood(node))
        subgraph.vs["label"]=subgraph.vs["name"]
        layout=subgraph.layout.kamada_kawai()
        igraph.plot(subgraph,layout=layout)

    def trajectory_graph(self,start_res=0,C_distance=5.4,non_C_distance=4.6):
        graph_trajectory=[]
        for ts in self.universe.trajectory:
            self.dagget_residue(start_res,C_distance,non_C_distance)
            graph_trajectory.append(self.graff)
        return graph_trajectory

    def analyze_trajectory_graph(self,trajectory_graph):
        for frame_num in range(len(trajectory_graph)):
            graph=trajectory_graph[frame_num]
            self.compute_graph_properties(graph)

    def graph_attribute_graphs(self,trajectory_graph):
        from scipy.interpolate import interp1d
        import matplotlib.pyplot as plt
        attributes=trajectory_graph[0].attributes()
        attributes_list={}
        for attribute in attributes:
            attributes_list[attribute]=[]
        for frame_num in range(len(trajectory_graph)):
            graph=trajectory_graph[frame_num]
            for attribute in attributes:
                attributes_list[attribute].append(graph[attribute])
        del attributes_list['laplacian']
        numplots=len(attributes_list.keys())
        fig=plt.subplot((numplots+1)/2,numplots/2)
        for key in attributes_list:
            plt.plot(attributes_list[key])
            plt.legend("{0}".format(key))
        plt.yscale('log')
        plt.show()
