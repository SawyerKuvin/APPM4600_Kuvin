##APPM 4600 - Homework 12
#Author: Sawyer Kuvin
#Date: 4/25/25

#-------Importing Libraries-------#
import numpy as np
from numpy import linalg as lg
import matplotlib.pyplot as plt
import scipy.special as sp
from scipy.linalg import lu_factor, lu_solve, inv

#-------Calculation Functions-------#
def hilbert(N):
    A = np.zeros((N,N))
    for i in range(N):
        for j in range(N):
            A[i,j] = 1/((i+1)+(j+1)-1)
    
    return(A)

def powerMethod(A,N,Nmax):
    q = np.random.rand(N,1)
    q = q/np.linalg.norm(q)
    lam = np.zeros((Nmax,1))

    for j in range(Nmax):
        z = np.matmul(A,q)
        q = z/np.linalg.norm(z)
        lam[j] = np.matmul(np.matmul(np.transpose(q),A),q)

    v = q
    return(lam,v)

#-------Problem Functions-------#
def Problem3():
    N = [4,8,12,16,20]
    print("Part a.")
    for n in N:
        A_mat = hilbert(n)
        #print(A_mat)

        Nmax = 100
        [lam, _] = powerMethod(A_mat,n,Nmax)

        j = 1
        lamMax = 0
        tol = 10**(-10)
        while abs(lam[j]-lam[j-1]) > tol:
            lamMax = lam[j]
            j = j+1

        print("n = ", n, "the largest eigen is ", lamMax, "and it took ", j, " iterations")

    print("Part b.")
    for n in N:  
        A_mat = hilbert(n)
        A_mat_inv = inv(A_mat)

        Nmax = 500
        [lam, _] = powerMethod(A_mat_inv,n,Nmax)
        #print(lam)

        j = 1
        lamMin = 0
        tol = 10**(-18)
        while abs(1/lam[j]-1/lam[j-1]) > tol:
            lamMin = 1/lam[j]
            j = j+1

        if N == 12 or N ==  16 or N == 20:
            lamMin = 1/lam[50]

        print("n = ", n, "the smallest eigen is ", lamMin, "and it took ", j, " iterations")



#-------Calling Problem Functions-------#
Problem3()