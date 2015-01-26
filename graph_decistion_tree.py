import numpy as np
import re
from scipy import stats
from sklearn import tree
from sklearn import ensemble

#mutation_list=['wt','f1174l','f1245c','f1245v','g1128a','i1170n','i1170s',
#               'i1171n','m1166r','r1192p','r1275q','t1151m','y1278s']
mutation_list={'wt':0,'d1349h':0,'r1231q':0,'f1174l':1,'r1275q':1,'f1245v':1}
directory='./'
start_res=1084
#c_dist_list=[4, 4.5 , 5, 5.5, 6]
#non_c_dist_list=[4, 4.5, 5, 5.5, 6]
c_dist_list=[4]
non_c_dist_list=[4]
dt=50

mutation_dataset={'mutation':[],'class':[],'feature_names':[],'data':[]}
for c_dist in c_dist_list:
    for non_c_dist in non_c_dist_list:
        for item in mutation_list.items():
            mutation,class_name = item
            out_file_name="{0}{1}_{2}_{3}_{4}".format(directory,mutation,c_dist,
                                                      non_c_dist,dt)
            fp=open(out_file_name+'.dat','r')
            lines=fp.readlines()
            fp.close()
            header=lines.pop(0) #remove header
            mutation_dataset['mutation'].append(mutation)
            mutation_dataset['class'].append(class_name)
            features=[]
            data_values=[]
            for line in range(len(lines)):
                if lines[line]=="\n" and line!=len(lines)-1:
                    split=re.split('\t',lines[line+1].rstrip())
                    features.append(split[0])
                    if split[1]!='nan':
                        data_values.append(float(split[1]))
                    else:
                        data_values.append(float(0))
            mutation_dataset['data'].append(data_values)
            mutation_dataset['feature_names']=features

def get_data():
    return mutation_dataset

X=np.array(mutation_dataset['data'])
y=np.array(mutation_dataset['class'])
forrest=ensemble.ExtraTreesClassifier(n_estimators=1000,random_state=0,
                                      criterion='entropy')
forrest.fit(X,y)
important=forrest.feature_importances_
std=np.std([tree.feature_importances_ for tree in forrest.estimators_],axis=0)
indices=np.argsort(important)[::-1]
for f in range(len(mutation_dataset['feature_names'])):
        print("{0:<4}{1} ({2})".format(f + 1, mutation_dataset['feature_names'][f], 
                                       important[indices[f]]))
