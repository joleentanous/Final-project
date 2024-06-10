import kmeansmodule as g
import pandas as pd
import sys
import math
import numpy as np

ITER = 300
error_messages = {
    "k": "Invalid number of clusters!",
    "iter": "Invalid maximum iteration!",
}


# needs to be fixed (like the kmeans in project1)
# errors might need to be handled by c errors system
# k, iter, eps, data1, data2
def main():
    k = iter = eps = data1 = data2 = None
    args = sys.argv
    k = args[1]
    
    if  len(sys.argv) >= 6:
        iter = args[2]
        eps = args[3]
        data1 = args[4]
        data2 = args[5]
        if (not iter.isdigit()): return (print(error_messages["iter"]))
        else: iter = int(iter)
        data = combine_data(data1, data2)
    else:
        iter = ITER
        eps = args[2]
        data1 = args[3]
        data2 = args[4]
        data = combine_data(data1, data2)
    
    eps= float(eps)
    #checks if the arguments are natural numbers
    if (not k.isdigit()): return (print(error_messages["k"]))
    k= int(k)
    #checks the validity of the arguments
    if (not 1<iter<1000): return print(error_messages["iter"])
    n = data.shape[0]
    d = data.shape[1]
    if (k<=1 or n<=k): return print(error_messages["k"])
    indices = (kmeansplus(k, n, d, iter, data))
    outputlist=(g.fit(k, n, d, iter, eps, indices, data.flatten().tolist()))
    output = arr_to_mtrx(outputlist,d,k)
    out = [[round(num, 4) for num in centroid] for centroid in output]
    list_str = ', '.join(map(str, indices))
    print(list_str)
    for row in out:
        print(",".join(map(str, row)))
    print("")





# joins two given vectors of length n/2 into one vector of length n
def combine_data(data1, data2):
    df1 = pd.read_csv(data1, header=None)
    df2 = pd.read_csv(data2, header=None)
    df1.columns=["key"]+df1.columns.tolist()[1:]
    df2.columns=["key"]+df2.columns.tolist()[1:]
    result_str=pd.merge(df1,df2,on="key")
    res = result_str.astype(float)
    res_sorted = res.sort_values(by="key")
    np_arr = res_sorted[res_sorted.columns[1:]].values
    return np_arr


# data is a 2D matrix of the dimensions nxd
def kmeansplus(k, n, d, iter, data):
    datacopy=data.copy()
    np.random.seed(0)
    ind = np.random.choice(n)
    rnd_clus = data[ind]
    clusters = np.zeros((k,d))
    clusters[0]=data[ind]
    indices = [ind] #represents the indices of the final chosen k clusters
    clus_len=1
    cur=0 #indicates the index of the last cluster added 
    for i in range(k-1):  # untill we get k clusters
        sum = 0
        DX_arr = []
        for row in data:
            DX = find_DX(clusters[:clus_len, :], row)
            sum+= DX
            DX_arr.append(DX)
        prob = [dx/sum for dx in DX_arr]
        new_clust_ind = np.random.choice(n, p=prob)
        indices.append(new_clust_ind)
        clusters[cur+1]=data[new_clust_ind]
        cur+=1
        clus_len+=1
    return indices


def euc_l2(v1, v2):
    dist = 0
    for i in range(len(v1)):
        dist += math.pow(v1[i] - v2[i], 2)
    return math.sqrt(dist)


def find_DX(clusters, v):
    DX = math.inf
    for row in clusters:
        distance = euc_l2(row, v)
        if distance < DX:
            DX = distance
    return DX

def arr_to_mtrx(outputlist,d,k):
    mtrx=[[0]*d for i in range(k)]
    for i in range(len(outputlist)):
        mtrx[i//d][i%d] = outputlist[i]
    return mtrx
            
    

#for checks
try:
    main()
except Exception as e:
    print("An Error Has Occured")
    print(e)
    
