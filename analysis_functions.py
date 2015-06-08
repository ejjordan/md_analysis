import numpy as np
import os
import gzip
import md_analysis as mda
from MDAnalysis.analysis.align import *

class hbonds_score(object):
    def __init__(self, mutation_list, hbonds_path=os.getcwd()):
        self.mutation_list=mutation_list
        self.hbonds_path=hbonds_path
        self.wt_act_donors=self.get_hbond_list(system_name='act_wt')
        self.wt_inact_donors=self.get_hbond_list(system_name='wt')

    def get_hbond_list(self,system_name,vmd_format=True):
        hbonds_file=open(self.hbonds_path + '/' + '{0}_hbonds.dat'.format(system_name),'rb')
        bonds=hbonds_file.readlines()
        hbonds_file.close()
        system_bonds=[]
        for i in range(len(bonds)):
            bonds[i]=bonds[i].strip(')\n')
            bonds[i]=bonds[i].strip('(')
            fields=bonds[i].split(', ')
            system_bonds.append(fields)

        if vmd_format==False:
            return system_bonds
            
        main_donor='N'
        main_acceptor='O'
        system_donors={}
        for i in range(len(system_bonds)):
            donor_idx=system_bonds[i][3]
            acceptor_idx=system_bonds[i][7]
            donor_id=system_bonds[i][2].strip("'")
            acceptor_id=system_bonds[i][6].strip("'")
            donor_heavy=system_bonds[i][4].strip("'")
            acceptor_heavy=system_bonds[i][8].strip("'")
            value=float(system_bonds[i][-1])
            if donor_heavy == main_donor:
                donor_loc='main'
            else:
                donor_loc='side'
            if acceptor_heavy == main_acceptor:
                acceptor_loc='main'
            else:
                acceptor_loc='side'
            if '{0}{1}-{2} {3}{4}-{5}'.format(donor_id,donor_idx,donor_loc,acceptor_id,
                                              acceptor_idx,acceptor_loc) in system_donors.keys():
                system_donors['{0}{1}-{2} {3}{4}-{5}'.format(donor_id,donor_idx,donor_loc,acceptor_id,
                                                               acceptor_idx,acceptor_loc)]+=value
            else:
                system_donors['{0}{1}-{2} {3}{4}-{5}'.format(donor_id,donor_idx,donor_loc,acceptor_id,
                                                               acceptor_idx,acceptor_loc)]=value
        return system_donors


    def find_baseline_altered_donors(self):
        all_donors=set(self.wt_inact_donors.keys()) | set(self.wt_act_donors.keys())
        altered_donors={}
        for key in all_donors:
            if key in self.wt_act_donors.keys():
                occupancy1=float(self.wt_act_donors[key])
            else:
                occupancy1=0
            if key in self.wt_inact_donors.keys():
                occupancy2=float(self.wt_inact_donors[key])
            else:
                occupancy2=0
            diff=abs(occupancy1-occupancy2)
            if diff>0:
                altered_donors[key]=diff
                #print "{0}\t{1:.4f}".format(key,diff)
        self.altered_donors=altered_donors
        self.bond_list=['ARG1275-side ASP1276-side','LYS1285-side ASP1160-side','ARG1279-side ASP1276-side',
                   'ARG1275-side ALA1274-main','ARG1279-side ASP1163-side','LYS1285-main GLY1304-main',
                   'ILE1277-main HSD1244-main','ARG1275-main ILE1246-main','ARG1275-side ASP1163-side',
                   'THR1310-side GLU1299-side','LYS1150-side ASP1270-side','SER1172-main ALA1168-main',
                   'LYS1285-side ASP1163-side','LEU1169-main LEU1165-main','TRP1295-main PRO1292-main',
                   'LYS1294-side ASP1249-side','ARG1373-side GLU1299-side','ARG1279-side GLN1159-side',
                   'ILE1171-main GLU1167-main','GLU1299-main GLU1299-side','TYR1283-main PHE1306-main',
                   'TYR1283-side MET1290-main','TYR1283-side GLY1287-main','ASN1175-main TYR1239-side',
                   'ALA1168-main PHE1164-main','ARG1284-side ASP1163-side','GLY1304-main GLU1299-main',
                   'TRP1295-side GLU1321-side']

    def print_baseline_altered_donors(self):
        for key in self.bond_list:
            if key in self.altered_donors.keys():
                print "found ", key, self.altered_donors[key]
            else:
                print "missing ", key

    def print_num_diff_hbonds(self,occupancy=False):
        self.find_baseline_altered_donors()
        for mutation in self.mutation_list:
            print mutation
            mut_inact_donors=self.get_hbond_list(system_name=mutation)
            count=0
            for key in self.bond_list:
                if key in self.wt_inact_donors.keys():
                     occupancy1=float(self.wt_inact_donors[key])
                else:
                    occupancy1=0
                if key in mut_inact_donors.keys():
                    occupancy2=float(mut_inact_donors[key])
                else:
                    occupancy2=0
                if occupancy:
                    print occupancy2*100
                    next
                delta_mut=abs(occupancy1 - occupancy2)
                if key in self.altered_donors.keys():
                    if self.altered_donors[key] > 0.4:
                        #print key, self.altered_donors[key]
                        delta_wt=float(self.altered_donors[key])
                        #print delta_wt, delta_mut, key
                    else:
                        next
                else:
                    next
                if delta_mut/delta_wt > 0.5:
                    #print key, delta_mut/delta_wt
                    count+=1
            if not occupancy:
                print count

class sasa_score(object):
    def __init__(self, mutation_list, hbonds_path=os.getcwd()):
        self.mutation_list=mutation_list
        self.probed_resids=[1096,1098,1170,1171,1174,1179,1239,1245,1271,1240]
        
    def read_sasa_file(self, system_name):
        variant_file=open('{0}_sasa.dat'.format(system_name),'rb')
        variant_sasa_lines=variant_file.readlines()
        variant_file.close()
        variant_sasa=[]
        variant_sasa_sum=0
        for line in variant_sasa_lines:
            data=line.split()
            variant_sasa.append(data)
            for res in self.probed_resids:
                if int(data[1]) == res:
                    variant_sasa_sum+=float(data[2])            
        return variant_sasa,variant_sasa_sum

    def print_sasa_scores(self,show_res=True,show_sums=False):
        wt_sasa_scores,wt_sasa_sum=self.read_sasa_file('wt')
        for mutant in self.mutation_list:
            print mutant
            sasa_scores,sasa_sum=self.read_sasa_file(mutant)
            for data in sasa_scores:
                for res in self.probed_resids:
                    if int(data[1]) == res:
                        if show_res:
                            print u"{0} {1}\t {2:1.3f}\u00B1{3:1.1f}".format(data[1],data[0],float(data[2]),
                                                                             float(data[3]))
                        else:
                            print u"{0:1.3f}\u00B1{1:1.1f}".format(float(data[2]),float(data[3]))
            if show_sums:
                print "delta SASA", (sasa_sum-wt_sasa_sum)
                print "total SASA", sasa_sum


class pca_score(object):

    def __init__(self, mutation_list, pca_path=os.getcwd()):
        self.mutation_list=mutation_list
        self.pca_path=pca_path
        self.wt_act_evals,self.wt_act_evecs=self.get_eigens_list(system_name='act_wt')
        self.wt_evals,self.wt_evecs=self.get_eigens_list(system_name='wt')
        self.num_top_pcs=50

    def get_eigens_list(self, system_name):
        system_file=gzip.open('{0}_eigenvecs.gz'.format(system_name),'rb')
        system_evecs=system_file.readlines()
        system_file.close()
        system_file=open('{0}_eigenvals.dat'.format(system_name),'rb')
        system_evals=system_file.readlines()
        system_file.close()
        for i in range(len(system_evals)):
            system_evals[i]=system_evals[i].strip('+0j)\n')
            system_evals[i]=float(system_evals[i].strip('('))
            system_evecs[i]=np.array(system_evecs[i].strip().split(),dtype=float)
        return system_evals,system_evecs

    def get_structure(self,directory,name):
        psf='{0}alk.dry.{1}.psf'.format(directory,name)
        dcd='{0}alk.40ns.dry.{1}.dcd'.format(directory,name)
        structure=mda.analyzer(psf,dcd,last_frame=0)
        return structure

    def difference_vector(self):
        wt_inact=self.get_structure("/home/joe/alk/analysis/40ns-dry-new-mutations/",'wt')
        wt_act=self.get_structure("/home/joe/alk/analysis/40ns-dry-new-mutations/",'act_wt')
        wt_inact_CA=wt_inact.universe.selectAtoms("name CA")
        wt_act_CA=wt_act.universe.selectAtoms("name CA")
        alignto(wt_act_CA,wt_inact_CA)
        diff_vector=(wt_act_CA.coordinates() - wt_inact_CA.coordinates()).flatten()
        #diff_vector=(wt_inact_CA.coordinates() - wt_act_CA.coordinates()).flatten()
        norm=np.linalg.norm(diff_vector)
        self.diff_vector=(diff_vector/norm)

    def make_pcs(self):
        self.difference_vector()
        for mutation in self.mutation_list:
            print "{0}".format(mutation)
            mutant_evals,mutant_evecs=self.get_eigens_list(mutation)
            overlap_sum=0
            cum_overlap=0
            cum_dot_diff=0
            for i in range(self.num_top_pcs):
                #cum_overlap+=(np.dot(self.wt_act_evecs[i],mutant_evecs[i]))
                cum_overlap+=(np.dot(self.diff_vector,mutant_evecs[i]))
                dot_diff=((np.dot(self.diff_vector,mutant_evecs[i])))*mutant_evals[i]
                cum_dot_diff+=dot_diff
                for j in range(self.num_top_pcs):
                    overlap_sum+=np.dot(self.wt_act_evecs[i],mutant_evecs[j])**2
            rmsip=np.sqrt(overlap_sum/self.num_top_pcs)
            self.princ_comp=mutant_evecs[0]
            dot_diff=(np.dot(self.diff_vector,self.princ_comp))*mutant_evals[0]
            #print "FO {0:7.4f}".format((np.dot(self.wt_act_evecs[0],mutant_evecs[0]))),
            print "FO {0:7.4f}".format((np.dot(self.diff_vector,mutant_evecs[0]))),
            print "CO {0:.4f}".format(cum_overlap), '\t',
            print "dot_diff {0:.4f}".format(dot_diff),
            print "dot_diff_sum {0:.4f}".format(cum_dot_diff),'\t',
            #print "rmsip {0:.4f}".format(rmsip),
            print "eval {0:.4f}".format(mutant_evals[0]),
            print "evsum {0:.4f}".format(np.sum(mutant_evals[0:self.num_top_pcs]))
