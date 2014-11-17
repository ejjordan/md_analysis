import sys
sys.path.insert(0,'/home/joe/git/md_analysis/')
import md_analysis as jma
import MDAnalysis as mda
import numpy as np
import pca as pca

mutation_list=['wt','f1174l','f1245c','f1245v','g1128a','i1170n','i1170s',
               'i1171n','m1166r','r1192p','r1275q','t1151m','y1278s']
#mutation_list=['r1192p','r1275q','t1151m','y1278s']
directory='./140ns-dry/'
#directory='/home/joe/alk/A1251T/'

for mutation in mutation_list:
    print "\nanalyzing {0} system".format(mutation)
    psf='{0}alk.dry.{1}.psf'.format(directory,mutation)
    dcd='{0}alk.140ns.aligned.dry.{1}.dcd'.format(directory,mutation)
#    psf="{0}alki.{1}.ions.psf".format(directory,mutation)
#    dcd="{0}alki.{1}.nvt-constrain.dcd".format(directory,mutation)
    out_file_name="{0}".format(mutation)
    PCA=pca.PCA(psf, dcd, outname=out_file_name)
    PCA.make_covariance()
    PCA.write_data()
    analysis=jma.analyzer(psf,dcd,out_file=out_file_name)
    analysis.get_rms()
    analysis.get_sasa()
    analysis.get_salts()
    analysis.get_hbonds()
    
