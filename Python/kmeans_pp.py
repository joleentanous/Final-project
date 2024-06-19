import symnmfmodule as g
import pandas as pd
import sys
import numpy as np
np.random.seed(0)

error_messages = {"k": "An Error Has Occurred"}
def read(data):
    df = pd.read_csv(data, header=None)
    np_array = df.to_numpy()
    return np_array


#initialization of the H matrix
def initH(W,k):
    n = W.shape[0]
    m = np.mean(W)
    H = np.random.uniform(0, 2 * np.sqrt(m / k), (n, k))
    return H

def symnmf(X, N, k):
    A = sym(X, N, k)
    D = ddg(A, N)
    W = norm(A, D, N)
    H = initH(W,k)
    return g.symnmf(W, H, N, k)

def printMatrix(matrix):
    for row in matrix:
        print(",".join(map(str, row)))
    print("")


def sym(X, N, d):
    return g.sym(X,N,d)


def ddg(A, N):
    return g.ddg(A, N)


def norm(A, D, N):
    return g.norm(A, D, N)




def main():
    k=goal=file_name = None
    args = sys.argv
    k= (args[1])
    goal= (args[2])
    file_name= (args[3])
    matrix=read(file_name)
    N=matrix.shape[0]
    if (k>=N or not k.isdigit()):
        return print(error_messages["k"])
    d=matrix.shape[1]
    if (goal=='symnmf'):
        res = symnmf(matrix,k)
        printMatrix(matrix)
    else:
        if (goal=='sym'):
            res = sym(matrix,N,d)
        elif (goal=='ddg'):
            A = sym(matrix, N, d)
            res = ddg(A,N)
        else:
            A = sym(matrix,N,d)
            D = ddg(A,N)
            res = norm(A, D, N)
        printMatrix(matrix)

