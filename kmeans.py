import math
import sys

EPSILON = 0.001

class Cluster():
    """
    A class representing a cluster of points in R^D.
    @type centroid: List[float]
    @param centroid: centroid of cluster.
    @type size: int
    @param size: size of cluster.
    @type cluster: List[List[float]]
    @param cluster: list of points in cluster.
    @type prev: List[float]
    @param prev: previous centroid.
    """
    def __init__(self, centroid):
        self.centroid = centroid
        self.size = 1
        self.cluster = [centroid]
        self.prev = [0 for i in range(len(self.centroid))]

    """
    Adds point to cluster list in cluster and updates size.
    @type point: List[float]
    @param point: point to add.
    """
    def update_cluster(self, point):
        self.cluster.append(point)
        self.size +=1

    """
    Removes all points from point list in cluster and updates size.
    """
    def clear_cluster(self):
        self.cluster = []
        self.size = 0

    """
    Calculates mean of all point coordinates and updates centroid, updated prev as well.
    @rtype: float
    @returns: the distance between previous and new centroid.
    """    
    def update_centroid(self):
        self.prev = self.centroid[:]
        n = len(self.prev)
        curr = [0 for i in range(n)]
        for i in range(n):
            for point in self.cluster:
                curr[i] += point[i]
            curr[i] /= self.size
        self.centroid = curr
        return self.calc_distance(self.prev)
    
    """
    Calculates euclidian disstance between cluster and a given point.
    @type point: List[float]
    @param point: point from data.
    @rtype: float
    @returns: the distance between previous and new centroid.
    """ 
    def calc_distance(self, point):
        dist = sum((point[i] - self.centroid[i]) ** 2 for i in range(len(point)))
        return math.sqrt(dist)


"""
Verifies argument is an int
@type arg: String
@param arg: A string of input from user to verify
@rtype: Touple(int, int)
@returns: First int return 1 if arg is not an integer and 1 if it is, second returns the arg in integer form.
""" 
def isInt(arg):
    try:
        return 0, int(arg)
    except:
        try:
            float_arg = float(arg)
            int_arg = int(arg[0 : arg.find(".")])
            if (float_arg) == int_arg:
                return 0, int_arg
            else:
                return 1, 0
        except:
            return 1, 0
    

"""
Parses file from txt file to a 2d matrix of coordinates.
@type filepath: String
@param filepath: path to txt file
@rtype: List[List[float]]
@returns: 2d matrix of all points from input.
"""
def get_points(filepath):
    with open(filepath, 'r') as file:
        data = [[float(num) for num in line.split(',')] for line in file]
    return data


"""
Verifies argument validity: correct number of arguments, correct input for K, iter, file.
Parses file from txt file to a 2d matrix of coordinates.
@type args: List[String]
@param args: arg[0] is file name, arg[1] is number of cluster k, arg[2] is either number of iterations or path to txt file, arg[3] is optional - if exists then its the filepath.
@rtype: List[int, int, List[List[float]]]
@returns: List[0] = k, List[1] = iter, List[2] = data OR [-1,0,0] if error.
""" 
def verify_data(args):
    if len(args) not in [3, 4]:
        print("An error has occurred!")
        return -1, 0, 0
    err_K,K = isInt(args[1])
    err_iter = 0
    if len(args) == 3:
        iterations = 200
    else:
        err_iter, iterations = isInt(args[2])
    path = args[2] if len(args) == 3 else args[3]

    if not path.endswith(".txt"):
        print("An error has occurred!")
        return -1, 0, 0
    
    data = get_points(path)

    if len(data) == 0:
        print("An error has occurred!")
        return -1, 0, 0
    if len(data) <= K or K <= 0 or err_K:
        print("Invalid number of clusters!")
        return -1, 0, 0
    if iterations >= 1000 or iterations <= 0 or err_iter:
        print("Invalid maximum iteration!")
        return -1, 0, 0
    
    return K, iterations, data


"""
Calculates distance between point and all cluster centroids and find the cluster with minimal distance.
Adds point to the point list of the closest cluster.
@type point: List[float]
@param point: point to calculate distance from.
@type cluster_list: List[Cluster]
@param cluster_list: list of clusters to find closest cluster from.
""" 
def find_closest_clus(point,cluster_list):
    min_dist = float('inf')
    index = 0
    K = len(cluster_list)
    for i in range(K):
        d = cluster_list[i].calc_distance(point)
        if d < min_dist:
            min_dist = d
            index = i
    cluster_list[index].update_cluster(point)


"""
Prints all centroid coordinates of clusters in cluster list.
@type cluster_list: List[Cluster]
@param cluster_list: list of all clusters.
""" 
def print_clusters(cluster_list):
    for cluster in cluster_list:
        target = ",".join([str(round(i, 4)) for i in cluster.centroid])
        print(target)


"""
Function is only relevant for the first iteration - it takes the first k lines in data and adds them to cluster list.
Afterwards it takes the rest of the lines in data and adds them to the closest cluster to the point.
@type K: int
@param k: number of clusters.
@type data: List[List[float]]
@param data: 2d matrix of all points from input.
@type cluster_list: List[Cluster]
@param cluster_list: cluster list
""" 
def initialize_clusters(K,data,cluster_list):
    for i, point in enumerate(data):
        if i < K:
            cluster_list.append(Cluster(point))
        else:
            find_closest_clus(point,cluster_list)


"""
Main logic of kmeans algorithm.
First initializes clusters and clears them.
Then runs a loop iterations times and calculates the centroids of the new clusters.
Loop also stops when the difference between all cluster centroids and prevs is less then EPSILON.
@type K: int
@param K: number of clusters.
@type iterations: int
@param iterations: number of iterations.
@type filepath: String
@param filepath: path to txt file.
@rtype: List[Cluster]
@returns: list of all clusters with their centroids.
"""
def kmeans(K, iterations, filepath):
    cluster_list = []
    data = get_points(filepath)
    initialize_clusters(K, data, cluster_list)
    for cluster in cluster_list:
        cluster.clear_cluster()

    while iterations > 0:
        for point in data:
            find_closest_clus(point,cluster_list)

        flag = True
        i = 0
        for cluster in cluster_list:
            i += 1
            if cluster.update_centroid() > EPSILON:
                flag = False            
            cluster.clear_cluster()

        if flag:
            break
        iterations -= 1

    return cluster_list


"""
Runs main logic of file and prints centroids of calculated k mean clusters.
@rtype: int
@returns: 0 if program runs successfully else 1.

""" 
def main():
    K, iterations, data = verify_data(sys.argv)
    if K == -1:
        return 1
    cluster_list = kmeans(K, iterations, data) 
    print_clusters(cluster_list)
    return 1


if __name__ == "__main__":
    main()
