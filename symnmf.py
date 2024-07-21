import symnmfmodule as g
import pandas as pd
import sys
import numpy as np
np.random.seed(0)

# symnmfmodule functions:
# All functions in the module accept a 1-d ==== PYTHON LIST ==== 
# Module functions:
# module_sym:  args: (X, N, k) where X is a py list (dim N*k), N, k integers
# module_ddg:  args: (X, N) where X is a py list (dim N*N), N integer
# module_norm:  args: (X, Y, N) where X and Y are  py lists (dim N*N), N integer
# module_symnmf:  args: (X, Y, N) where X and Y are  py lists (X with dim N*N, Y with dim N*k), N, k integers


def error():
    print("An Error Has Occurred, python")
    exit()

def read(data):
    print(2)
    df = pd.read_csv(data, header=None)
    print(3)
    np_array = df.to_numpy()
    print(np_array)
    return np_array


#initialization of the H matrix
def initH(N,W,k):
    m = np.mean(W)
    H = np.random.uniform(0, 2 * np.sqrt(m / k), (N, k))
    H = H.flatten().tolist()
    return H

def symnmf(X, N, k):
    A = sym(X, N, k)
    D = ddg(A, N)
    W = norm(A, D, N)
    H = initH(N,W,k)
    return g.module_symnmf(W, H, N, k)


def printMatrix(matrix , N, K):
    matrix = np.array(matrix).reshape(N, K)
    matrix = [[round(num, 4) for num in row] for row in matrix]
    for row in matrix:
        print(",".join(map(str, row)))
    print("")



def sym(X, N, d):
    return g.module_sym(X,N,d)


def ddg(A, N):
    return g.module_ddg(A, N)


def norm(A, D, N):
    return g.module_norm(A, D, N)




def main():
    goals = ['symnmf', 'sym', 'ddg', 'norm']
    k = goal = file_name = None
    args = sys.argv
    if (len(args) != 4):
        error()
    k = (args[1])
    goal = (args[2])
    file_name = (args[3])
    print(1)
    matrix = read(file_name)
    N = matrix.shape[0]
    try:
        k = int(k)
        goal = str(goal)
        file_name = str(file_name)
    except:
        error()

    if (k >= N) or (k <= 1) or (goal not in goals):
        error()



    d = matrix.shape[1]
    matrix = matrix.flatten().tolist()
    print(matrix)
    print(4)
    if goal == 'symnmf':
        res = symnmf(matrix,N,k)
        print(res)
        printMatrix(res, N, k)
    else:
        if goal == 'sym':
            res = sym(matrix,N,d)
        elif goal == 'ddg':
            A = sym(matrix, N, d)
            res = ddg(A,N)
        elif goal == 'norm':
            A = sym(matrix,N,d)
            D = ddg(A,N)
            res = norm(A, D, N)
        else:
            error()
        printMatrix(res, N, N)

# try:
#     main()
# except Exception as e:
#     error()

main()
