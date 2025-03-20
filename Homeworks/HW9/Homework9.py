##APPM 4600 - Homework 9
#Author: Sawyer Kuvin
#Date: 3/21/25


#Importing Commonly Used Libraries
import numpy as np
from numpy import linalg as lg
import matplotlib.pyplot as plt
import scipy.special as sp
from scipy.linalg import lu_factor, lu_solve, inv

#Problem 2 - Computation
def Problem2():
    A_ls = np.array([[581,105],[105,454]])
    b_ls = np.array([421,247])
    #print(inv(A_ls))
    x_vec = np.matmul(inv(A_ls) , b_ls)

    print("The x-vector is ", x_vec)


#Calling Problem Functions
Problem2()