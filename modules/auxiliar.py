# --------------------------------------------------------------------------------------------------
# Title: Grav Codes
# Author: Edson Alonso Falla Luza
# Description: Source codes 
# Collaboratores: Rodrigo Bijani
# --------------------------------------------------------------------------------------------------

# -------- Import Python internal libraries ---------
import math

def sec(angle):
    sec = (1/math.cos(angle))
    return sec

def cossec(angle):
    cossec = (1/math.sin(angle))
    return cossec

def cotg(angle):
    cotg = math.cos(angle)/math.sin(angle)
    return cotg