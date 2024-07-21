import sys
import numpy as np
import pandas as pd
from sklearn.metrics import silhouette_score
import kmeans
import symnmf

ITER = 200


def error():
    print("An Error Has Occurred, analysis")
    exit()


def read(data):
    df = pd.read_csv(data, header=None)
    np_array = df.to_numpy()
    return np_array



def symnmf_clustering(matrix, N, k):
    print("symnmf_clus")
    D_matrix = matrix.flatten().tolist()
    print(1)
    symnmf_H = np.array(symnmf.symnmf(D_matrix, N, k))
    print(2)
    indices = symnmf_H.argmax(axis=1)
    print(indices)
    return indices

def kmeans_clustering(matrix, N, k, d):
    print("kmeans clus")
    centroids = kmeans.k_means(k, N, d, ITER, matrix)   
    indices = []
    print(centroids)
    for datapoint in matrix:
        indices.append(kmeans.find_closest_centroid_index(centroids, datapoint))
    print(indices)
    return indices


if __name__ == "__main__":
    print("hi")
    args = sys.argv
    if (len(args) != 3):
        error()
    k = (args[1])
    file_name = (args[2])
    try:
        k = int(k)
        file_name = str(file_name)
    except:
        error()

    matrix = read(file_name)
    N = matrix.shape[0]
    d = matrix.shape[1]



    kmeans_clusters = kmeans_clustering(matrix, N, k ,d)  # calculating the centroids using kmeans
    symnmf_clusters = symnmf_clustering(matrix, N ,k)  # calculating the centroids using symnmf

    print("nmf: %.4f" % silhouette_score(matrix, symnmf_clusters))
    print("kmeans: %.4f" % silhouette_score(matrix, kmeans_clusters))





# import sys
# import numpy as np
# from sklearn.metrics import silhouette_score
# from sklearn.cluster import KMeans
# from symnmf import symnmf, sym, ddg, norm  # Assuming you have these functions in symnmf.py

# def load_data(file_name):
#     """Load data from a text file."""
#     return np.loadtxt(file_name, delimiter=',')

# def compute_symnmf(file_name, k):
#     """Compute SymNMF and return the silhouette score."""
#     # Load data
#     X = load_data(file_name)
    
#     # Perform SymNMF
#     H = symnmf(X, k)  # Ensure this function returns the matrix H
    
#     # Compute similarity matrix
#     W = norm(X)
    
#     # Compute cluster assignments from H
#     cluster_assignments = np.argmax(H, axis=1)
    
#     # Calculate silhouette score
#     score = silhouette_score(X, cluster_assignments)
#     return score

# def compute_kmeans(file_name, k):
#     """Compute K-means and return the silhouette score."""
#     # Load data
#     X = load_data(file_name)
    
#     # Perform K-means
#     kmeans = KMeans(n_clusters=k, random_state=0).fit(X)
    
#     # Calculate silhouette score
#     score = silhouette_score(X, kmeans.labels_)
#     return score

# def main():
#     if len(sys.argv) != 3:
#         print("Usage: python3 analysis.py <k> <file_name>")
#         sys.exit(1)
    
#     k = int(sys.argv[1])
#     file_name = sys.argv[2]
    
#     # Compute silhouette scores
#     nmf_score = compute_symnmf(file_name, k)
#     kmeans_score = compute_kmeans(file_name, k)
    
#     # Print results
#     print(f"nmf: {nmf_score:.4f}")
#     print(f"kmeans: {kmeans_score:.4f}")

# if __name__ == "__main__":
#     main()
