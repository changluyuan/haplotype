#coding=UTF-8
import itertools 
import sys,os
#from numpy.linalg import inv,qr
import numpy as np
import gzip
import re
from collections import OrderedDict
if (len(sys.argv)<>4):
   print("Usage: "+ sys.argv[0]+ '  sample_readbed.gz  CGMAP haplotype_site.cfg')
   exit(-1) 

infile1 = gzip.open(sys.argv[1])
infile2 = open(sys.argv[2])
file1name = os.path.basename(sys.argv[1])
file2name = os.path.basename(sys.argv[3])
name1 = file1name.strip().split('_readbed')[0]
name2 = file2name.strip().split('.cfg')[0]
np.seterr(divide='ignore',invalid='ignore')
np.set_printoptions(threshold='nan')
if re.search(r'minus',sys.argv[2]):
    flag = 'minus'
if re.search(r'plus',sys.argv[2]):
    flag = 'plus'
outfile = open(name1+'_'+name2+'_'+flag+'.txt','w')
haplotype_site = OrderedDict()
####20180424 change haplotype_site[site] to haplotype_site[hap_site]
with open(sys.argv[3]) as cfg:
    for line in cfg:
        list1 = line.strip().split('\t')
        haplotype = list(map(int,list1[0].strip().split(',')))
        site = tuple(list1[1].strip().split(',')[1:])
        chrom = list1[1].strip().split(',')[0]
        if len(haplotype) == len(site):
            hap_site = tuple(zip(site,haplotype))
            haplotype_site[hap_site] = chrom
        else:
            print('Error haplotype: ' + line.strip())

##############################################################################
##hap is a list,containing haplotypes. e.g. [1,1],[-1,-1]
def vector_generator(hap,site):
    vector = []
    baseline = []
    for i in range(0,length):
        baseline.append(0)
    continu = len(hap)
    for elem in range(0,continu):
        if (site[0],site[elem+1]) not in CGMAP_rev:
            sys.stderr.write(','.join(site)+'\t'*2 +str(site[elem+1])+' is out of CGMAP region'+'\n')
            baseline = []
            break
        relative_pos = CGMAP_rev[(site[0],site[elem+1])]
        baseline[relative_pos] = hap[elem]
    vector = tuple(baseline)
    vector_array = np.array(vector)
    #vector_array_T = vector_array.T
    return(vector_array)  

def generate_new_array(array):
    complement_zero = np.zeros((array.shape[0],length - array.shape[1]))
    array_new = np.hstack((array,complement_zero))
    return(array_new)

def compare_to_V(hapsite):
## scan every region(about 300bp)
    site_tmp = list(zip(*list(hapsite))[0])
    chro = haplotype_site[hapsite]
    site = [chro] + site_tmp
    hap = list(zip(*hapsite)[1])
    continu =len(hap)
    data_array_pre = np.array([])
    group_num = len(Group)
    i = 0
    for key in Group:
        i += 1
        valid_length = CGMAP_range[(key[0],key[1],key[2])]
        data_array = np.array(Group[key],dtype=np.float64)
        if site[0] == key[0] and int(site[1]) - int(key[1]) >= 600 and int(site[1]) - int(key[1]) <= 1200:
            data_array_left = data_array[:,valid_length:]
            data_array_pre = generate_new_array(data_array_left)
        elif site[0] == key[0] and int(site[1]) - int(key[1]) >= 0 and int(key[2]) - int(site[1]) >=0:
            data_array_valid = data_array[:,:valid_length]
            data_array_pro = generate_new_array(data_array_valid)
            if data_array_pre.size:
                data_array_all = np.vstack((data_array_pre,data_array_pro))
                data_array_pre = np.array([])
            else:
                data_array_all = data_array_pro
            vector_array_T = vector_generator(hap,site) 
            if not vector_array_T.size:
                print(','.join(map(str,hap))+'\t'+','.join(map(str,site))+'  do not have vector')
                break
            result = np.dot(data_array_all,vector_array_T)
            result_abs = np.dot(abs(data_array_all),abs(vector_array_T))
## scan every vector_T column(every haplotype)
            right = len(np.argwhere(result==continu))
            valid = len(np.argwhere(result_abs==continu))
            all_num =  len(np.argwhere(result_abs != 0))
#####calculate the ratio of specific haplotype in all the reads containing all haplotype region
            ratio = float(right)/valid if valid != 0 else 0
#####calculate the ratio of specific haplotype in all the reads containing all or half haplotype region
            ratio_all = float(right)/all_num if all_num != 0 else 0
            outfile.write(site[0]+'\t'+':'.join(site[1:])+'\t'+':'.join(map(str,hap))+'\t'+str(all_num)+'\t'+str(valid)+'\t'+str(right)+'\t'+str(ratio)+'\n')  
            break
        if i == len(Group):
            print(','.join(map(str,hap))+'\t'+','.join(map(str,site))+'  not in group(readbed.gz)')
###################################################################################
def main():
##CGMAP
    global CGMAP,CGMAP_range,length,CGMAP_rev,Group
    CGMAP = OrderedDict()
    CGMAP_rev = OrderedDict()
    CGMAP_range =OrderedDict()
    length = 0
    for line in infile2:
        list1 = line.strip().split('\t')
        chrom = list1[0]
        region_start = list1[1]
        region_end = list1[2]
        absolute_pos = list1[3]
        relative_pos = int(list1[4])
        CGMAP_range[(chrom,region_start,region_end)] = relative_pos + 1
        CGMAP[(chrom,region_start,region_end,relative_pos)] = absolute_pos
        CGMAP_rev[(chrom,absolute_pos)] = relative_pos
      
###################################################################################
## read the inputfile, save epiread to a dictionary
    Group = OrderedDict()
    for line in infile1:
        list1 = line.strip().split('\t')
        if list1[0] == 'chr':
            continue
        else:
            row = list1[3].split(',')
            length = len(row)
            chrom = list1[0]
            start = int(list1[1])
            end = int(list1[2])
            if(list1[0],list1[1],list1[2]) not in Group:
                Group[(list1[0],list1[1],list1[2])] =[row]
            else:
                Group[(list1[0],list1[1],list1[2])].append(row)
    for site in haplotype_site:
        compare_to_V(site)
if __name__ =="__main__":
    main()

