# --------------------------------------------------------------------------------------------------
# Title: Grav Codes
# Author: Edson Alonso Falla Luza
# Description: Source codes 
# Collaboratores: Rodrigo Bijani
# --------------------------------------------------------------------------------------------------

# -------- Import Python internal libraries ---------
import numpy as np
import math
from math import sqrt
import networkx as nx
from scipy.spatial.distance import mahalanobis

# -------- L1 Norm (sum norm) ---------

def l1dist(x_coord,y_coord):
    '''
    This function takes two vectors, x_coord and y_coord, and returns a matrix were the element in the ij position is the distance (considering the sum norm, ou l1 norm) betwen the point (x_coord[i],y_coord[i]) and (x_coord[j],y_coord[j]).
    
    Inputs:
    x_coord - numpy array 
    y_coord - numpy array
    
    Output:
    dl1 - numpy array - Matrix of distances
    '''
    
    #Stablishing the error condition
    tamx = np.shape(x_coord)[0]           
    tamy = np.shape(x_coord)[0]
    if tamx != tamy:
        raise ValueError("All inputs must have same length!")
        
    #Calculating and savingn the distances
    #tam = ( math.factorial(tamx) )/( 2*(math.factorial(tamx-2)) ) #2 choises over 'tamx'(or 'tamy') posibilities 
    distl1_matrix = np.zeros((tamx,tamy))
                                    
    for i in range(tamx):
        for j in range(tamx):
            distl1_matrix[i][j] = abs(x_coord[i] - x_coord[j]) + abs(y_coord[i] - y_coord[j])                                      
            
    return distl1_matrix
    

# -------- L2 Norm (euclidians norm) ---------

def l2dist(x_coord,y_coord):
    '''
    This function takes two vectors, x_coord and y_coord, and returns a matrix were the element in the ij position is the distance (considering the euclidian norm, ou l2 norm) betwen the point (x_coord[i],y_coord[i]) and (x_coord[j],y_coord[j]).
    
    Inputs:
    x_coord - numpy array 
    y_coord - numpy array
    
    Output:
    dl1 - numpy array - Matrix of distances
    '''
    
    #Stablishing the error condition
    tamx = np.shape(x_coord)[0]           
    tamy = np.shape(x_coord)[0]
    if tamx != tamy:
        raise ValueError("All inputs must have same length!")
        
    #Calculating and savingn the distances
    #tam = ( math.factorial(tamx) )/( 2*(math.factorial(tamx-2)) ) #2 choises over 'tamx'(or 'tamy') posibilities
    distl2_matrix = np.zeros((tamx,tamy))
                                    
    for i in range(tamx):
        for j in range(tamx):
            distl2_matrix[i][j] = ( (x_coord[i] - x_coord[j])**2 + (y_coord[i] - y_coord[j])**2 )**(1/2)                                     
            
    return distl2_matrix


# -------- Graph L2 Distance Function ---------

def distgraph(G):
    '''
    This function recives a graph with the edges already weighed and returns a value (float) of the square of the sum of the weight of all knots to the average of the weights (phi) as the average itself (dmG).
    
    Inputs:
    G = graph already weighed
    
    Output:
    
    phi = float
    dmG = float
    '''
    # get the weights of the undirected Graph:
    dm1 = []
    for (u, v, wt) in G.edges.data('weight'): 
        dm1.append(wt)
   
    # computing the mean of the MST:
    dmG = np.mean( np.array(dm1) )

    # compute the variance of the distances of the MST:
    phi = 0.0
    for i in range( len(dm1) ):
        phi += (dm1[i] - dmG)**2
    
    return phi,dmG

# -------- Graph L1 Distance Function ---------

def distgraphl1(G):
    '''
    This function recives a graph with the edges already weighed and returns a value (float) of the sum of the weight of all knots to the average of the weights (phi) as the average itself (dmG).
    
    Inputs:
    G = graph already weighed
    
    Output:
    
    phi = float
    dmG = float
    '''
    # get the weights of the undirected Graph:
    dm1 = []
    for (u, v, wt) in G.edges.data('weight'): 
        dm1.append(wt)
   
    # computing the mean of the MST:
    dmG = np.mean( np.array(dm1) )

    # compute the variance of the distances of the MST:
    phi = 0.0
    for i in range( len(dm1) ):
        phi += abs(dm1[i] - dmG)
    
    return phi,dmG

# -------- Graph Mahalanobis Distance Function ---------

def distgraphmaha(x,y):
    '''
    This function recives two vectors, meaning the x coordinates and the y coordinates, computes a graph weighted with the distances betwen the knots(x,y), calculates the minimum spaning tree and returns the mahalanobis distances of the vertices in the MST.
    Inputs:
    x - array like (n,1)
    y - arrat like (n,1)
    Output:
    
    m - float - mahalanobis distance
    '''
    #creates the graph and the MST:
    S = nx.Graph()
    for i in range(len(x)):
        S.add_node(i ,pos=(x[i],y[i]))
    for j in range(len(x)):
        S.add_edge(i,j,weight=l2dist(x,y)[i][j])
    TS = nx.minimum_spanning_tree(S)
    
    # get the weights of the undirected Graph:
    dm1 = []
    u   = []
    v   = []
    for (i, j, wt) in TS.edges.data('weight'): 
        u.append(x[i])
        v.append(y[j])
    N = np.array( (u,v) )
    C = np.cov( N.T )
    invC = np.linalg.inv(C)
    m = mahalanobis(u,v,invC)
   
    
    return m