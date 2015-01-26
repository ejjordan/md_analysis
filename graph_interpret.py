import numpy as np
import re
import scipy.stats as stats

#mutation_list=['wt','f1174l','f1245c','f1245v','g1128a','i1170n','i1170s',
#               'i1171n','m1166r','r1192p','r1275q','t1151m','y1278s']
mutation_list=['wt','f1174l','r1275q','d1349h','r1231q']
directory='./'
start_res=1084
#c_dist_list=[4, 4.5 , 5, 5.5, 6]
#non_c_dist_list=[4, 4.5, 5, 5.5, 6]
c_dist_list=[4]
non_c_dist_list=[4]
dt=50

mutation_info={}
for c_dist in c_dist_list:
    for non_c_dist in non_c_dist_list:
        for mutation in mutation_list:
            #print "\nanalyzing {0}_{1}_{2}_{3} system".format(mutation,c_dist,
            #                                                non_c_dist,dt)
            out_file_name="{0}{1}_{2}_{3}_{4}".format(directory,mutation,c_dist,
                                                      non_c_dist,dt)
            fp=open(out_file_name+'.dat','r')
            lines=fp.readlines()
            fp.close()
            header=lines.pop(0) #remove header
            feature_list=[]
            for line in range(len(lines)):
                if lines[line]=="\n" and line!=len(lines)-1:
                    split=re.split('\t',lines[line+1].rstrip())
                    almost_list=lines[line+2].rstrip(']\n')
                    almost_list=almost_list.lstrip('[')
                    value_list=re.split(', ',almost_list)
                    value_list=[float(x) for x in value_list]
                    split.append(np.array(value_list))
                    feature_list.append(split)
            mutation_info[mutation]=feature_list


p_vals={}
control=mutation_info.pop('wt')
for key in mutation_info.keys():
    #print "\n{0}".format(key)
    p_vals[key]=[]
    value=mutation_info[key]
    for i in range(len(value)):
        feature_name=value[i][0]
        t,p=stats.ttest_ind(control[i][3],value[i][3],axis=0)
        #print "{0}:\nt stat = {1} p val = {2}".format(feature_name,t,p)
        p_vals[key].append([feature_name,t,p])

values=p_vals.values()
for i in range(len(values)):
    for j in range(len(values)/2):
        if i!=j:
            print 
            for k in range(len(values[i])):
                print values[i][k][0],values[i][k][2],values[j][k][2]
