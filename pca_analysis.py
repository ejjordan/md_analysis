import numpy as np
import gzip

mutation_list=['f1174l','f1245v','f1245c','i1170n','i1170s','i1171n','y1278s',
               'r1192p','m1166r','r1275q','t1151m','g1128a','g1286r','a1200v',
               'd1349h','r1231q','d1270g']
#mutation_list=['v1229m'] #fix this one

wt_file=gzip.open('wt_peter_eigenvecs.gz','rb')
wt_evecs=wt_file.readlines()
wt_file.close()
wt_file=open('wt_peter_eigenvals.dat','rb')
wt_evals=wt_file.readlines()
wt_file.close()

num_top_pcs=9
for mutation in mutation_list:
    print mutation,
    top_wt_pcs=[]
    top_mutant_pcs=[]
    mutant_file=gzip.open('{0}_peter_eigenvecs.gz'.format(mutation),'rb')
    mutant_evecs=mutant_file.readlines()
    mutant_file.close()
    mutant_file=open('{0}_peter_eigenvals.dat'.format(mutation),'rb')
    mutant_evals=mutant_file.readlines()
    mutant_file.close()
    for i in range(0,num_top_pcs):
        wt_pc=np.array(wt_evecs[i].strip().split(),dtype=float)
        mutant_pc=np.array(mutant_evecs[i].strip().split(),dtype=float)
        if i==1:
            print "first overlap {0:7.4f}".format(np.dot(wt_pc,mutant_pc)),
        top_wt_pcs.append(wt_pc)
        top_mutant_pcs.append(mutant_pc)

    overlap_sum=0
    cum_overlap=0
    for i in range(num_top_pcs):
        cum_overlap+=np.abs(np.dot(top_wt_pcs[i],top_mutant_pcs[i]))
        for j in range(num_top_pcs):
            overlap_sum+=np.dot(top_wt_pcs[i],top_mutant_pcs[j])**2
    rmsip=np.sqrt(overlap_sum/10)
    print "rmsip {0:.4f}".format(rmsip),
    print "cumulative overlap {0:.4f}".format(cum_overlap),
    print "eval {0}".format(mutant_evals[0]),
