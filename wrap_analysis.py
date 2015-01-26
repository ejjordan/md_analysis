#import sys
#sys.path.insert(0,'/home/joe/git/md_analysis/')
import md_analysis as jma
import MDAnalysis as mda
import numpy as np
import pca as pca

#mutation_list=['wt','f1174l','f1245c','f1245v','g1128a','i1170n','i1170s',
#               'i1171n','m1166r','r1192p','r1275q','t1151m','y1278s']
mutation_list=['a1200v','d1270g','d1349h','g1286r','r1231q']
#mutation_list=['v1229m'] #fix this one
directory='../40ns-dry/'
num_frames=2000 #want 40 ns worth of trajectory

for mutation in mutation_list:
    print "\nanalyzing {0} system".format(mutation)
    psf='{0}alk.dry.{1}.psf'.format(directory,mutation)
    dcd='{0}alk.40ns.dry.{1}.dcd'.format(directory,mutation)
#    psf="{0}alki.{1}.ions.psf".format(directory,mutation)
#    dcd="{0}alki.{1}.nvt-constrain.dcd".format(directory,mutation)
    out_file_name="{0}_peter".format(mutation)
    PCA=pca.PCA(psf, dcd, align_selection="name CA and not ((resid 1125:1129) or (resid 1157:1174) or (resid 1271:1292))", pca_selection="name CA", outname=out_file_name, last_frame=num_frames)
    PCA.make_covariance()
    PCA.write_data()
"""
    analysis=jma.analyzer(psf,dcd,out_file=out_file_name,last_frame=num_frames)
    analysis.get_rms()
    analysis.get_sasa()
    analysis.get_salts()
    analysis.get_hbonds()
"""
