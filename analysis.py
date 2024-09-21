import numpy as np
from sklearn.metrics import silhouette_score
import math
import sys
import kmeans
import symnmf


"""
Converts command line arguments to usable data
@rtype: Tuple(int, string)
@returns: A tuple of the parsed data containing K and filepath to data
"""
def convert_args():
    if len(sys.argv) != 3:
        print("An Error has occurred!")
        sys.exit(1)
    Err_k, k = kmeans.isInt(sys.argv[1])
    if Err_k:
        print("An Error has occurred!")
        sys.exit(1)
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
            dist = kmeans.calc_distance(point,centroids[i])
            if dist < min_dist:
                min_dist = dist
                index = i
        res.append(index)
    return res


"""
Runs the symnmf algorithm on a given points set
@type K: int
@param K: The number of clusters
@type points: List[List[float]]
@param points: A list of points
@rtype: List[List[float]]
@returns: A list of H matrix
"""
def do_symnmf(K,points):
    res = symnmf.symnmf(K, points)
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
    np_points = np.array(points)

    kmeans_centroids = do_kmeans(k,300,filepath)
    kmeans_res = get_points_centroid_index(points,kmeans_centroids)

    symnmf_H = do_symnmf(np_points, k)
    symnmf_res = index_of_max_in_lists(symnmf_H)

    print("nmf: %.4f" % silhouette_score(points, symnmf_res))
    print("kmeans: %.4f" % silhouette_score(points, kmeans_res))


if __name__ == "__main__":
    main()
