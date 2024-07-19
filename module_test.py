import symnmfmodule as g
import sys
import math

def printMatrix1DAs2D(matrix, N, d):
    for i in range(N):
        for j in range(d):
            print(f"{matrix[i * d + j]:.6f}", end=" ")
        print()
    print()

ITER = 300
def main():
    N = 3
    d = 2

    X = [
        1.0, 2.0,
        3.0, 4.0,
        5.0, 6.0
    ]

    similarityMatrix = g.module_sym(X, N, d)

    print("Similarity Matrix:")
    printMatrix1DAs2D(similarityMatrix, N, N)

    
    ddg_matrix = g.module_ddg(similarityMatrix, N)
    print("DDG Matrix:")
    printMatrix1DAs2D(ddg_matrix, N, N)
    # printMatrix(ddg_matrix, N, N)

    print("Normalized Matrix:")
    normalized_matrix = g.module_norm(similarityMatrix, ddg_matrix, N)
    printMatrix1DAs2D(normalized_matrix, N, N)


    N = 3
    k = 2
    W = [
        0.000000, 0.707105, 0.000006,
        0.707105, 0.000000, 0.707105,
        0.000006, 0.707105, 0.000000
    ]

    H = [
        0.5, 0.6,
        0.3, 0.8,
        0.9, 0.1
    ]
    optimizedH = g.module_symnmf(W, H, N, k)
    print("SYMNMF Matrix:")
    printMatrix1DAs2D(optimizedH, N, k)


main()

def printMatrix(matrix, rows, cols):
    for i in range(rows):
        for j in range(cols):
            print(f"{matrix[i * cols + j]:.2f}", end=" ")
        print()
    print()


