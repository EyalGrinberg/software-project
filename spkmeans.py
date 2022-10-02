import spkmeansmodule
import sys
from sys import argv
import numpy as np
import pandas as pd
import math

def kmeansPP(k, max_iter, eps, t_mat): # k++ algorithm

    N = len(t_mat) #rows
    dims = len(t_mat[0]) #columns

    #k-means++
    data = t_mat
    sum_dists = 0
    i = 0
    np.random.seed(0)
    new_cent = np.random.randint(0, N) # choose first centroid randomly
    dists_probs = np.full((2,N) , math.inf) # create distances and probabilities array
    centroids = np.zeros((1,dims)) + data[new_cent] # create centroids array and insert first centroid 
    indices = [0 for i in range(k)] 
    indices[0] = new_cent # first indice in result list
    while (i < k-1): # main loop
        dists_probs[0] = np.minimum(dists_probs[0] , np.sum(np.power((data - centroids[i]) , 2) , axis = 1)) # keep min distance
        sum_dists = np.sum(dists_probs[0] , axis = 0) # sum of distances from calculated centroid
        np.true_divide(dists_probs[0] , sum_dists , out = dists_probs[1]) # point probability = divide distance from centroid by sum of distances
        new_cent = np.random.choice(N , p = dists_probs[1]) # choose new index according to the probabilities
        i += 1
        indices[i] = (int(new_cent)) # add new index
        centroids = np.vstack([centroids , data[new_cent]]) # add new centroid
    return indices

def normalize_u(u_mat): 
    t_mat = u_mat
    for i in range(len(u_mat)):
        sum_j = 0
        for j in range(len(u_mat[0])):
            sum_j += (u_mat[i][j] ** 2)
        for j in range(len(u_mat[0])):
            if sum_j == 0:
                t_mat[i][j] = 0
            else:
                t_mat[i][j] = u_mat[i][j] / (sum_j ** 0.5)
    return t_mat

def final_output(k, output_mat_c, goal): #final adjustments and printing
    if goal == 5: # only for "spk" - retrieve k centroids and print k indices
        init_centroids = kmeansPP(k , 300 , 0 , output_mat_c) # get initial k++ centroids
        for ind in init_centroids: # cast to float for C api
            ind = float(ind)   
        centP = [str(a) for a in init_centroids] # for printing
        print(','.join(centP)) #print indices from k++
        output_mat_c = output_mat_c.tolist() # convert np.array to list for C-api
        output_mat_c = spkmeansmodule.kmeansCapi(output_mat_c , init_centroids , k) #get final centroids from kmeans C   
    # print output matrix by requested format   
    new_res = [[0 for i in range(len(output_mat_c[0]))] for j in range(k)] 
    for i in range(k): 
        for j in range(len(output_mat_c[0])):
            new_res[i][j] = output_mat_c[i][j] 
            if j == len(output_mat_c[0]) - 1:
                print(format(new_res[i][j], ".4f"))
            else:
                print(format(new_res[i][j], ".4f"), end = '')
                print("," , end = '')
    print("")  

def extract_from_file(file_name): # retrieve data points 
    data_points = pd.read_csv(file_name, sep=",", header=None)
    N = len(data_points)
    dims = len(data_points.columns)
    res = np.array(data_points)
    res = res.tolist()
    return res

#checking arguments
if len(argv) != 4: #check num of args 
    print("Invalid Input!")
    sys.exit(1)
try:
    k = int(argv[1])
    if k < 0:
        print("Invalid Input!")
        sys.exit(1)
except:
    print("Invalid Input!")
    sys.exit(1)
goal = argv[2]
if (goal!="spk") and (goal!="wam") and (goal!="ddg") and (goal!="lnorm") and (goal!="jacobi") :
    print("Invalid Input!")
    sys.exit(1)

file_name = argv[3]

data_points = extract_from_file(file_name)

if goal == "spk":
    goal = 5 # encoding goal to int
    c_V_mat = spkmeansmodule.spkCapi(data_points) #receive jacobi's output
    k_sortedV = spkmeansmodule.eigenGapCapi(c_V_mat) #receive k from heuristic and sorted V from jacobi [k, sortedV]
    if k == 0: # obtain k from eigenGap heuristic, else k remains as given in input
        k = k_sortedV[0]
    eigen_vectors = k_sortedV[1]
    eigen_vectors = eigen_vectors[1:] #leave only the vectors
    if k >= len(eigen_vectors): # k >= N
        print("Invalid Input!")
        sys.exit(1)
    else:
        u_mat = np.array(eigen_vectors)
        u_mat = np.delete(u_mat, np.s_[k:], 1) # take only first k eigen vectors
        t_mat = normalize_u(u_mat)
        final_output(k, t_mat, goal) # send T to next stages of the algorithm: k++ and kmeans
else: #goal != spk 
    # encoding goal to int
    if goal == "wam": goal = 1
    if goal == "ddg": goal = 2
    if goal == "lnorm": goal = 3
    if goal == "jacobi": goal = 4
    c_output = spkmeansmodule.pyToCtoPy(data_points, goal) #main func will return wam/ddg/lnorm/jacobi output
    k = len(c_output) # adjust k for printing dimensions
    final_output(k, c_output, goal) # send for printing
