__author__ = 'slhuang'
import os
from subprocess import Popen, PIPE
from itertools import tee, izip


def rocchio_nosampling(up_exp,dn_exp):
	matrixF = '/workspace/code/May28/matrix/feature_gene_scale'
	#'/workspace/project1/LINCS/networks/feature_gene_scale'
	up_exp = '/workspace/code/Dec11/exp/'+up_exp
	dn_exp='/workspace/code/Dec11/exp/'+dn_exp
	prog = '/workspace/code/May28/rocchio_modular/rocchio'
	proc = Popen([prog,'-matrixF',matrixF,'-expF',up_exp,'-expF2',dn_exp,'-weighted', '-sortG', '-earlyT'],stdout=PIPE)
	proc.wait()
	
'''
def rocchio_nosampling(up_exp):
	matrixF = '/workspace/knowdata/KnowNets/KnowNet_0.3/LINCS/MCF7_data_sets/MCF7_drug.level4.mat.txt'
	up_exp = '/workspace/knowdata/KnowNets/KnowNet_0.3/LINCS/MCF7_data_sets/drug_exp_sets/'+up_exp
	
	prog = '/workspace/code/May28/rocchio_modular/rocchio'
	proc = Popen([prog,'-matrixF',matrixF,'-expF',up_exp, '-delimiter', 'TRANSPOSE','-weighted', '-sortG', '-earlyT'])
	proc.wait()
'''

def rocchio_sampling(up_exp,dn_exp):
	matrixF = '/workspace/code/May28/matrix/feature_gene_scale'
	#'/workspace/project1/LINCS/networks/feature_gene_scale'
	up_exp = '/workspace/code/Dec11/exp/'+up_exp
	dn_exp='/workspace/code/Dec11/exp/'+dn_exp
	prog = '/workspace/code/May28/rocchio_modular/rocchio'
	proc = Popen([prog,'-matrixF',matrixF,'-expF',up_exp,'-expF2',dn_exp,'-weighted', '-sampl'])
	proc.wait()

'''
def rocchio_sampling(up_exp):
	matrixF = '/workspace/knowdata/KnowNets/KnowNet_0.3/LINCS/MCF7_data_sets/MCF7_drug.level4.mat.txt'
	up_exp = '/workspace/knowdata/KnowNets/KnowNet_0.3/LINCS/MCF7_data_sets/drug_exp_sets/'+up_exp
	
	prog = '/workspace/code/May28/rocchio_modular/rocchio'
	proc = Popen([prog,'-matrixF',matrixF,'-expF',up_exp, '-delimiter', 'TRANSPOSE','-weighted', '-sampl'])
	proc.wait()
'''
def rocchio_sampling_opt(up_exp,dn_exp):
	matrixF = '/workspace/code/May28/matrix/feature_gene_scale'
	#'/workspace/project1/LINCS/networks/feature_gene_scale'
	up_exp = '/workspace/code/Dec11/exp/'+up_exp
	dn_exp='/workspace/code/Dec11/exp/'+dn_exp
	
	prog = '/workspace/code/May28/rocchio_modular/rocchio'
	proc = Popen([prog,'-matrixF',matrixF,'-expF',up_exp,'-expF2',dn_exp,'-weighted', '-samplOpt'])
	proc.wait()

'''
def rocchio_sampling_opt(up_exp):
	matrixF = '/workspace/knowdata/KnowNets/KnowNet_0.3/LINCS/MCF7_data_sets/MCF7_drug.level4.mat.txt'
	up_exp = '/workspace/knowdata/KnowNets/KnowNet_0.3/LINCS/MCF7_data_sets/drug_exp_sets/'+up_exp
	
	prog = '/workspace/code/May28/rocchio_modular/rocchio'
	proc = Popen([prog,'-matrixF',matrixF,'-expF',up_exp, '-delimiter', 'TRANSPOSE','-weighted', '-samplOpt'])
	proc.wait()
'''
def rocchio_baseline(up_exp,dn_exp):
	matrixF = '/workspace/code/May28/matrix/feature_gene_scale'
	#'/workspace/project1/LINCS/networks/feature_gene_scale'
	up_exp = '/workspace/code/Dec11/exp/'+up_exp
	dn_exp='/workspace/code/Dec11/exp/'+dn_exp
	prog = '/workspace/code/May28/rocchio_modular/rocchio'
	proc = Popen([prog,'-matrixF',matrixF,'-expF',up_exp,'-expF2',dn_exp,'-weighted']) #,'-Fconsider',str(100),'-topK',str(100),'-vertical'])
	proc.wait()
'''	
def rocchio_baseline(up_exp):
	matrixF = '/workspace/knowdata/KnowNets/KnowNet_0.3/LINCS/MCF7_data_sets/MCF7_drug.level4.mat.txt'
	up_exp = '/workspace/knowdata/KnowNets/KnowNet_0.3/LINCS/MCF7_data_sets/drug_exp_sets/'+up_exp
	
	prog = '/workspace/code/May28/rocchio_modular/rocchio'
	proc = Popen([prog,'-matrixF',matrixF,'-expF',up_exp, '-delimiter', 'TRANSPOSE','-weighted'])
	proc.wait()
'''
def rocchio_notransform(up_exp,dn_exp):
	matrixF = '/workspace/code/May28/matrix/feature_gene_scale'
	#'/workspace/project1/LINCS/networks/feature_gene_scale'
	up_exp = '/workspace/code/Dec11/exp/'+up_exp
	dn_exp='/workspace/code/Dec11/exp/'+dn_exp
	prog = '/workspace/code/May28/rocchio_modular/rocchio'
	proc = Popen([prog,'-matrixF',matrixF,'-expF',up_exp,'-expF2',dn_exp,'-weighted','-notransform'])
	proc.wait()
	
# ./rocchio -matrixF /workspace/code/May28/matrix/cmb_combat_lm_corrected_exprData_Jan15_2016 -expF o -notransform -delimiter COMMA
	
def rocchio_sampling_opt_sortF_TOP1000(up_exp,dn_exp):
	matrixF = '/workspace/code/May28/matrix/feature_gene_scale'
	#'/workspace/project1/LINCS/networks/feature_gene_scale'
	up_exp = '/workspace/code/Dec11/exp/'+up_exp
	dn_exp='/workspace/code/Dec11/exp/'+dn_exp
	prog = '/workspace/code/May28/rocchio_modular/rocchio'
	proc = Popen([prog,'-matrixF',matrixF,'-expF',up_exp,'-expF2',dn_exp,'-weighted', '-samplOpt','-sortF','-Fconsider',str(20000000)])
	proc.wait()

def rocchio_sampling_opt_sortF_vertical_TOP1000(up_exp,dn_exp):
	matrixF = '/workspace/code/May28/matrix/feature_gene_scale'
	#'/workspace/project1/LINCS/networks/feature_gene_scale'
	up_exp = '/workspace/code/Dec11/exp/'+up_exp
	dn_exp='/workspace/code/Dec11/exp/'+dn_exp
	prog = '/workspace/code/May28/rocchio_modular/rocchio'
	proc = Popen([prog,'-matrixF',matrixF,'-expF',up_exp,'-expF2',dn_exp,'-weighted', '-samplOpt','-vertical','-sortF','-Fconsider',str(20000000)])
	proc.wait()
'''
def rocchio_sampling_opt_sortF_TOP1000(up_exp):
	matrixF = '/workspace/knowdata/KnowNets/KnowNet_0.3/LINCS/MCF7_data_sets/MCF7_drug.level4.mat.txt'
	up_exp = '/workspace/knowdata/KnowNets/KnowNet_0.3/LINCS/MCF7_data_sets/drug_exp_sets/'+up_exp
	
	prog = '/workspace/code/May28/rocchio_modular/rocchio'
	proc = Popen([prog,'-matrixF',matrixF,'-expF',up_exp, '-delimiter', 'TRANSPOSE','-weighted', '-samplOpt','-sortF','-Fconsider',str(1000)])
	proc.wait()
'''
def pairwise(iterable):
    "s -> (s0,s1), (s1,s2), (s2, s3), ..."
    a, b = tee(iterable)
    next(b, None)
    return izip(a, b)
    
def rocchio_run():
    dirs=sorted(os.listdir('/workspace/code/Dec11/exp'))
    #dirs=sorted(os.listdir('/workspace/knowdata/KnowNets/KnowNet_0.3/LINCS/MCF7_data_sets/drug_exp_sets'))
    print dirs
    for i in range(len(dirs)/2):
    	#if i >= 10:
    	#print i, dirs[i] 
    	#rocchio_sampling_opt_sortF_TOP1000(dirs[i])
    		#rocchio_nosampling(dirs[i])
    		#rocchio_baseline(dirs[i])
    	#rocchio_sampling(dirs[i])
    	#rocchio_sampling_opt(dirs[i])
    	print i,dirs[2*i+1],dirs[2*i]
    #for up_exp,dn_exp in pairwise(dirs):
        #if dir.endswith('.stp'):
            #continue
        #print up_exp,dn_exp
        #pre(dir)
        #rocchio_nosampling(dirs[2*i+1],dirs[2*i])
        #rocchio_sampling(dirs[2*i+1],dirs[2*i])
        #rocchio_sampling_opt(dirs[2*i+1],dirs[2*i])
        #rocchio_baseline(dirs[2*i+1],dirs[2*i])
        #rocchio_sampling_opt_sortF_TOP1000(dirs[2*i+1],dirs[2*i])
        #rocchio_baseline(dirs[2*i+1],dirs[2*i])
        rocchio_sampling_opt_sortF_vertical_TOP1000(dirs[2*i+1],dirs[2*i])
        #rocchio_notransform(dirs[2*i+1],dirs[2*i])
        #baseImp(dir,3)
        #baseImp(dir,2)



def run():
    rocchio_run()
               
if __name__ == '__main__':
    run()
    