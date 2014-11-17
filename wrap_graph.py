import sys
sys.path.insert(0,'/home/joe/git/md_analysis/')
import graphs as gp
import MDAnalysis as mda
import numpy as np

#mutation_list=['wt','f1174l','f1245c','f1245v','g1128a','i1170n','i1170s']
#mutation_list=['i1171n','m1166r','r1192p','r1275q','t1151m','y1278s']
mutation_list=['a1200v','r1231q','g1286r','v1229m']
directory='./40ns-dry/'
start_res=1084
c_dist_list=[4, 4.5, 5, 5.5, 6]
non_c_dist_list=[4, 4.5, 5, 5.5, 6]
dt=50
is_psf=True

for c_dist in c_dist_list:
    for non_c_dist in non_c_dist_list:
        for mutation in mutation_list:
            print "\nanalyzing {0}_{1}_{2}_{3} system".format(mutation,c_dist,
                                                            non_c_dist,dt)
            psf='alk.dry.{0}.psf'.format(mutation)
            dcd='alk.40ns.dry.{0}.dcd'.format(mutation)
            out_file_name="{0}_{1}_{2}_{3}".format(mutation,c_dist,non_c_dist,dt)
            universe=mda.Universe(directory+psf,directory+dcd)
            graf=gp.graphs(universe)
            trajectory_graph=graf.trajectory_graph(start_res=start_res,C_distance=
                                                   c_dist,non_C_distance=non_c_dist,
                                                  skip=dt,psf=is_psf)
            graf.analyze_trajectory_graph(trajectory_graph)
            graph_dict=graf.graph_attribute_graphs(trajectory_graph,display=False,
                                                   write_file=out_file_name+'.png')
            fp=open(out_file_name+'.dat','w')
            fp.write("start_res={0}\tC_distance={1}\tnon_C_distance={2}\tskip={3}\n\n"
                     .format(start_res,c_dist,non_c_dist,dt))
            for key in graph_dict:
                mean=np.mean(graph_dict[key])
                std=np.std(graph_dict[key])
                fp.write("{0}\t{1}\t{2}\n{3}\n\n".format(key,mean,std,graph_dict[key]))
            fp.close()
