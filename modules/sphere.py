# --------------------------------------------------------------------------------------------------
# Title: Grav Codes
# Author: Edson Alonso Falla Luza
# Description: Source codes 
# Collaboratores: Rodrigo Bijani
# --------------------------------------------------------------------------------------------------

# -------- Import Python internal libraries ---------
import numpy as np

def gz2D(x, z, sphere):
    '''    
    This function calculates the vertical component of gravity attraction produced by a solid point source. This is a Python 
    implementation for the subroutine presented in Blakely (1995). On this function, there are 
    received the value of the initial and final observation points (X and Y) and the properties 
    of the sphere.
       
    Inputs:
    sphere - numpy array - elements of the sphere
    sphere[0,1,2]
    sphere = list sphere[x_center(meters), z_center(meters), mass(kg)]
    
    Output:
    gz - numpy array - vertical component for the gravity in mGal. Size of gz is the same as x and z observations    
    '''
    
    # Stablishing some conditions
    if x.shape != z.shape:
        raise ValueError("All inputs must have same shape!")
    
    # Setting the initial value
    gz = 0.
    
    # Setting position-vector: 
    dx = sphere[0] - x
    dz = sphere[1] - z
    
    # properties of the sphere:
    mass = sphere[2]
    
    # Definition for some constants
    G = 6.673e-11 # SI
    si2mGal = 100000.0 
    
    # Compute the distance
    r = np.sqrt(dx**2 + dz**2)
    
    # Compute the vertical component 
    gz = mass * dz / (r**3)
    gz *= G*si2mGal
    
    # Return the final outpu
    return gz
###########################################################################################################################

def volfont(raio,rho):
    '''
    Nao esquecer de comentar o que essa funcao faz!
    '''
    V = (4.0/3.0) * np.pi * raio**3 * rho
    massa = V * rho
    return raio,rho,V,massa