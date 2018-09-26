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
            distl2_matrix[i][j] = sqrt( (x_coord[i] - x_coord[j])**2 + (y_coord[i] - y_coord[j])**2 )                                     
            
    return distl2_matrix


# -------- Graph Distance Function ---------

def distgraph(G):
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
        phi += (dm1[i] - dmG)**2
    
    return phi,dmG