##APPM 4600 - Lab 8 (In Lab Code Attempt)
#Author: Sawyer Kuvin
#Date: 3/4/25
#Questions attempted/completed: 
# (~ partially completed/attempted) 

import numpy as np
import time
from numpy import linalg as lg
from numpy.linalg import inv 
from numpy.linalg import norm
import matplotlib.pyplot as plt
import scipy.special as sp
from scipy.linalg import lu_factor, lu_solve

#-------Given/Reused Function-------#
def LineConstruct(x_0,x_1,f_x,alpha): #Pre-Lab
    if abs(x_0-x_1)<1e-10:
        print("Cannot define line")
    m = (f_x(x_1)-f_x(x_0))/(x_1-x_0)
    line = lambda x: f_x(x_0) + m*(x-x_0)
    f_alpha = line(alpha)
    return(f_alpha)

#-------Pre-Lab/Lab Questions-------#
def Prelab():
    f_x  = lambda x: x+2
    x_0 = 1
    x_1 = 23
    alpha = 4
    f_alpha = LineConstruct(x_0,x_1,f_x,alpha)
    print(f_alpha)



#-------Calling Problem Functions-------#
#Pre-Lab Test
Prelab()



