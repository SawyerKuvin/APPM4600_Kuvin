##APPM 4600 - Lab 8 (In Lab Code Attempt)
#Author: Sawyer Kuvin
#Date: 3/4/25
#Questions attempted/completed: 2, 3.1, 3.2, ~3.3
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

def  eval_lin_spline(xeval,Neval,a,b,f,Nint):
    '''create the intervals for piecewise approximations'''
    xint = np.linspace(a,b,Nint+1)
    '''create vector to store the evaluation of the linear splines'''
    yeval = np.zeros(Neval)
    for j in range(Nint):
        '''find indices of xeval in interval (xint(jint),xint(jint+1))'''
        '''let ind denote the indices in the intervals'''
        atmp = xint[j]
        btmp= xint[j+1]
        # find indices of values of xeval in the interval
        ind= np.where((xeval >= atmp) & (xeval <= btmp))
        xloc = xeval[ind]
        n = len(xloc)
        '''temporarily store your info for creating a line in the interval of
        interest'''
        fa = f(atmp)
        fb = f(btmp)
        yloc = np.zeros(len(xloc))
        for kk in range(n):
        #use your line evaluator to evaluate the spline at each location
            yloc[kk] = LineConstruct(atmp,btmp,f,xloc[kk])
        yeval[ind] = yloc
    return yeval

def create_natural_spline(yint,xint,N):

#    create the right  hand side for the linear system
    b = np.zeros(N+1)
#  vector values
    h = np.zeros(N+1)
    for i in range(1,N):
        hi = xint[i]-xint[i-1]
        hip = xint[i+1] - xint[i]
        b[i] = (yint[i+1]-yint[i])/hip - (yint[i]-yint[i-1])/hi
        h[i-1] = hi
        h[i] = hip

#  create the matrix A so you can solve for the M values
    A = np.zeros((N+1,N+1))
    A[0][0] = 1.0
    for j in range(1,N):
        A[j][j-1] = h[j-1]/6
        A[j][j] = (h[j]+h[j-1])/3
        A[j][j+1] = h[j]/6
        A[N][N] = 1
    Ainv = inv(A)
    M = Ainv.dot(b)
    
#  Create the linear coefficients
    C = np.zeros(N)
    D = np.zeros(N)
    for j in range(N):
        C[j] = yint[j]/h[j]-h[j]*M[j]/6
        D[j] = yint[j+1]/h[j]-h[j]*M[j+1]/6
    return(M,C,D)
       
def eval_local_spline(xeval,xi,xip,Mi,Mip,C,D):
# Evaluates the local spline as defined in class
# xip = x_{i+1}; xi = x_i
# Mip = M_{i+1}; Mi = M_i

    hi = xip-xi
   
    yeval = 0
    return yeval 
    
    
def  eval_cubic_spline(xeval,Neval,xint,Nint,M,C,D):
    
    yeval = np.zeros(Neval+1)
    
    for j in range(Nint):
        '''find indices of xeval in interval (xint(jint),xint(jint+1))'''
        '''let ind denote the indices in the intervals'''
        atmp = xint[j]
        btmp= xint[j+1]
        
#   find indices of values of xeval in the interval
        ind= np.where((xeval >= atmp) & (xeval <= btmp))
        xloc = xeval[ind]

# evaluate the spline
        yloc = eval_local_spline(xloc,atmp,btmp,M[j],M[j+1],C[j],D[j])
#   copy into yeval
        yeval[ind] = yloc

    return(yeval)
    

#-------Pre-Lab/Lab Questions-------#
def Prelab():
    f_x  = lambda x: x+2
    x_0 = 1
    x_1 = 23
    alpha = 4
    f_alpha = LineConstruct(x_0,x_1,f_x,alpha)
    print(f_alpha)

def Problem3_1():
    f = lambda x: np.exp(x)
    a = 0
    b = 1
    ''' create points you want to evaluate at'''
    Neval = 100
    xeval = np.linspace(a,b,Neval)
    ''' number of intervals'''
    Nint = 10
    '''evaluate the linear spline'''
    yeval = eval_lin_spline(xeval,Neval,a,b,f,Nint)
    ''' evaluate f at the evaluation points'''
    fex = f(xeval)
    plt.figure()
    plt.plot(xeval,fex,'ro-')
    plt.plot(xeval,yeval,'bs-')
    plt.legend(["Original Function", "Linear Spline"])
    plt.title("Estimation Comparison")
    plt.show()
    err = abs(yeval-fex)
    plt.figure()
    plt.plot(xeval,err,'ro-')
    plt.title("Linear Spline Error")
    plt.show()

def Problem3_2():
    f = lambda x: 1/(1+(10*x)**2)
    a = -1
    b = 1
    ''' create points you want to evaluate at'''
    Neval = 100
    xeval = np.linspace(a,b,Neval)
    ''' number of intervals'''
    Nint = 10
    '''evaluate the linear spline'''
    yeval = eval_lin_spline(xeval,Neval,a,b,f,Nint)
    ''' evaluate f at the evaluation points'''
    fex = f(xeval)
    plt.figure()
    plt.plot(xeval,fex,'ro-')
    plt.plot(xeval,yeval,'bs-')
    plt.legend(["Original Function", "Linear Spline"])
    plt.title("Estimation Comparison")
    plt.show()
    err = abs(yeval-fex)
    plt.figure()
    plt.plot(xeval,err,'ro-')
    plt.title("Linear Spline Error")
    plt.show()

def Problem3_3():
    f = lambda x: np.exp(x)
    a = 0
    b = 1
    
    ''' number of intervals'''
    Nint = 3
    xint = np.linspace(a,b,Nint+1)
    yint = f(xint)

    ''' create points you want to evaluate at'''
    Neval = 100
    xeval =  np.linspace(xint[0],xint[Nint],Neval+1)

#   Create the coefficients for the natural spline    
    (M,C,D) = create_natural_spline(yint,xint,Nint)

#  evaluate the cubic spline     
    yeval = eval_cubic_spline(xeval,Neval,xint,Nint,M,C,D)
    
    
    ''' evaluate f at the evaluation points'''
    fex = f(xeval)
        
    nerr = norm(fex-yeval)
    print('nerr = ', nerr)
    
    plt.figure()    
    plt.plot(xeval,fex,'ro-',label='exact function')
    plt.plot(xeval,yeval,'bs--',label='natural spline') 
    plt.legend
    plt.show()
     
    err = abs(yeval-fex)
    plt.figure() 
    plt.semilogy(xeval,err,'ro--',label='absolute error')
    plt.legend()
    plt.show()



#-------Calling Problem Functions-------#
#Pre-Lab Test
Prelab()

#Problem 3.1
Problem3_1()

#Problem 3.2
Problem3_2()

#Problem 3.3
Problem3_3()


