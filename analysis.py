import numpy as np
from sklearn.metrics import silhouette_score
import math
import sys
import kmeans
import symnmfmodule
from tools import isInt, error

np.random.seed(1234)

"""
Converts command line arguments to usable data
@rtype: Tuple(int, string)
@returns: A tuple of the parsed data containing K and filepath to data
"""
def convert_args():
    if len(sys.argv) != 3:
        error()
    k = isInt(sys.argv[1])
    filepath = sys.argv[2]
    return k, filepath


"""
Runs the kmeans algorithm on a given file
@type k: int
@param k: The number of clusters
@type iterations: int
@param iterations: The number of iterations to run the algorithm
@type filepath: string
@param filepath: The path to the file containing the data
@rtype: List[List[float]]
@returns: A list of cluster centroids

"""
def do_kmeans(k,iterations,filepath):
    res = kmeans.kmeans(k,iterations,filepath)
    centroids = [cluster.centroid for cluster in res]
    return centroids


"""Calculates euclidian disstance between cluster and a given point.
    @type pointA: List[float]
    @param pointA: point from data.
    @type pointB: List[float]
    @param pointB: point from data.
    @rtype: float
    @returns: the distance between 2 points
""" 
def distance(pointA, pointB):
    dist = sum((pointA[i] - pointB[i]) ** 2 for i in range(len(pointA)))
    return math.sqrt(dist)


"""
Finds the centroid index for each point
@type points: List[List[float]]
@param points: A list of points
@type centroids: List[List[float]]
@param centroids: A list of centroids
@rtype: List[int]
@returns: A list of centroid indices
"""
def get_points_centroid_index(points,centroids):
    res = []
    for point in points:
        min_dist = float('inf')
        index = 0
        for i in range(len(centroids)):
            dist = distance(point,centroids[i])
            if dist < min_dist:
                min_dist = dist
                index = i
        res.append(index)
    return res


"""
Runs the symnmf algorithm on a given points set
@type K: int
@param K: The number of clusters
@type points: NUMPY ARRAY
@param points: A list of points
@rtype: List[List[float]]
@returns: A list of H matrix
"""
def do_symnmf(K, points):
    N = len(points)
    W = symnmfmodule.norm(points)
    W_np = np.array(W)
    initial_H_np = np.random.uniform(low=0, high=(2*np.sqrt(np.average(W_np)/K)), size=(N, K))
    initial_H = initial_H_np.tolist()
    res = 0
    try:
        res = symnmfmodule.symnmf(initial_H, W)
    except:
        error()
    return res


"""
Finds the index of the maximum value in each list
@type H: List[List[float]]
@param H: A list of lists
@rtype: List[int]
@returns: A list of indices
"""
def index_of_max_in_lists(H):
    res = []
    for i in range(len(H)):
        if H[i] == []:
            res.append(None)
            continue
        max_val = -math.inf
        index = 0
        for j in range(len(H[i])):
            if H[i][j] > max_val:
                max_val = H[i][j]
                index = j
        res.append(index)
    return res


"""
Main function to compare SymNMF and KMeans using silhouette score.    
This function reads the input arguments, loads the data, runs both KMeans and SymNMF,
and then computes and prints the silhouette score for each clustering algorithm.
The silhouette score measures the quality of clustering by comparing within-cluster
distance against between-cluster distance.
"""
def main():
    k, filepath = convert_args()
    points = kmeans.get_points(filepath)
    kmeans_centroids = do_kmeans(k,300,filepath)
    kmeans_res = get_points_centroid_index(points,kmeans_centroids)

    symnmf_H = do_symnmf(k, points)
    symnmf_res = index_of_max_in_lists(symnmf_H)

    print("nmf: %.4f" % silhouette_score(points, symnmf_res))
    print("kmeans: %.4f" % silhouette_score(points, kmeans_res))

"""
This is the entry point of the program.
"""
if __name__ == "__main__":
    main()
