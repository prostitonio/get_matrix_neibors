import struct
from progress.bar import Bar
import argparse
from collections import Counter
import matplotlib.pyplot as plt
import numpy as np
import gc 
from matplotlib import cm
from matplotlib.colors import LightSource
import sys
from scipy.optimize import curve_fit
#name = r"voro.vor"

def get_frame(fobj,num):
    all_vol = []
    all_nei = []
    nei     = []
    for i in range(num):
        vol, = struct.unpack("f",fobj.read(4))
        num_nei, = struct.unpack("i",fobj.read(4))
        for j in range(num_nei):
            one_nei, = struct.unpack("i",fobj.read(4))
            nei.append(one_nei)
        nei = np.array(nei)
        all_nei.append(nei)
        all_vol.append(vol)
        nei = []
    return np.array(all_vol,dtype=object), np.array(all_nei,dtype=object)

def get_count_frame(name):
    try:
        fobj = open(name,'rb')
        num, = struct.unpack_from("i",fobj.read(4))
        num_nei,=struct.unpack_from("i",fobj.read(4+4),offset=4)
        count = 0
        while True:
            try:
                num_nei, = struct.unpack_from("i",fobj.read(8+num_nei*4),offset=4*(num_nei+1))
                count +=1
            except:
                return int( (count+1)/num )
    except:
        print("can not open *.vor file")
        return 0

def init(name):
    try:
        f = open(name);f = [i for i in f]
        print("open file ",name)
        f = np.array([[int(i[:5]),i[5:8]] for i in f[2:-1]]).T # array [ [num mol] , [name mol] ] 
        temp = Counter(f[1]) #[ [name mol , count mol] , ...]
        #print(temp)


        name_mol_first_type   = f[1][0]
        name_mol_second_type  = f[1][-1]
        count_atom_first_type = int(temp[name_mol_first_type])
        count_mol_first_type  = int(f[0][count_atom_first_type-1])
        count_all_mol         = int(f[0][-1])
        count_mol_second_type = count_all_mol - count_mol_first_type 


        print("system consists of ", name_mol_first_type," and ",name_mol_second_type)
        print("count all mol",str(count_all_mol) )
        print("count ",name_mol_first_type," mol: ", str(count_mol_first_type) )
        print("count ",name_mol_second_type," mol: ",str(count_mol_second_type) )
        print("concentration: " , name_mol_first_type , " x=%.3f" % (count_mol_first_type/count_all_mol)     )
        print("concentration: " , name_mol_second_type , " x=%.3f" % (count_mol_second_type/count_all_mol)     )

        #f.close()

        return count_mol_first_type


    except:
        print("can not open *.gro file ",name)
        sys.exit()
        print("test")



parser = argparse.ArgumentParser(description='Process some integers.')
parser.add_argument('-g', dest='gro', action='store', default="conf.gro", help='input *.gro file (default: conf.gro)')
parser.add_argument('-v', dest='vor', action='store', default="voro.vor", help='input *.vor (default: voro.vor)')
parser.add_argument('--verbose', dest='verbose', action='store_true', default=False)
args = parser.parse_args()

count_mol_first_type =            init(args.gro)
count_frame          = get_count_frame(args.vor)

#'''
M_max = 50
N_M_max = 50 

bar_nei_A = np.zeros([M_max,N_M_max])
bar_vol_A = np.zeros([M_max,N_M_max])
bar_nei_B = np.zeros([M_max,N_M_max])
bar_vol_B = np.zeros([M_max,N_M_max])

l = count_mol_first_type

with open(args.vor, 'rb') as fobj:
    num, = struct.unpack("i",fobj.read(4))
    all_vol_0 = []; all_vol_1 = [];
    all_nei_0 = []; all_nei_1 = [];

    bar = Bar('Processing', max=count_frame)
    for i in range(count_frame):
        vol , nei = get_frame(fobj,num)
        vol_0, vol_1 = vol[:l], vol[l:]
        nei_0, nei_1 = nei[:l], nei[l:]
        #if(i%10 == 0):print("THF_"+str(round(l/num,3))+":" +str(i))

        a_0=[];a_1=[]
        for i , j in enumerate(nei_0):
            m,n_m = j[j<l].shape[0]  ,  j[j>=l].shape[0] 
            bar_nei_A[m,n_m] +=1 
            bar_vol_A[m,n_m] +=vol_0[i] 
        for i , j in enumerate(nei_1):
            m,n_m = j[j<l].shape[0]  ,  j[j>=l].shape[0] 
            bar_nei_B[m,n_m] +=1 
            bar_vol_B[m,n_m] +=vol_1[i] 
        
        nei_0=[];nei_1=[]
        gc.collect()
        bar.next()
    bar.finish()


np.savetxt("solution_%.3f_nei_THF" % round(l/num,3) ,bar_nei_A)
print("save solution_%.3f_nei_THF" % round(l/num,3))
np.savetxt("solution_%.3f_vol_THF" % round(l/num,3),bar_vol_A)
print("save solution_%.3f_vol_THF" % round(l/num,3))
np.savetxt("solution_%.3f_nei_H2O" % round(l/num,3),bar_nei_B)
print("save solution_%.3f_nei_H2O" % round(l/num,3))
np.savetxt("solution_%.3f_vol_H2O" % round(l/num,3),bar_vol_B)
print("save solution_%.3f_vol_H2O" % round(l/num,3))
#'''
