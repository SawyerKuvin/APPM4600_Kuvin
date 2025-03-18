##APPM 4600 - Lab 10 (In Lab Code Attempt)
#Author: Sawyer Kuvin
#Date: 3/18/25
#Questions attempted/completed: 2, 3.2, ~3.3
# (~ partially completed/attempted) 


import matplotlib.pyplot as plt
import numpy as np
import numpy.linalg as la
import math
from scipy.integrate import quad

#--------Built Functions-------#
def eval_legendre(N,xi):
    p_vec = np.zeros((N+1))
    if N == 0:
        phi0 = 1
        p_vec[0] =phi0
        return(p_vec)
    else:
        phi0 = 1
        phi1 = xi
        p_vec[0] = phi0
        p_vec[1] = phi1
        
        for n in range(1,N):
            phin_plus_1 = (1/(n+1)) * (((2*n)+1)*xi*p_vec[n]-n*p_vec[n-1])
            p_vec[n+1] = phin_plus_1

        return(p_vec)
    
def eval_last_legendre(N,xi):
    p = eval_legendre(N,xi)
    return(p[N])

def eval_legendre_expansion(f,a,b,w,n,x): 

#   This subroutine evaluates the Legendre expansion

#  Evaluate all the Legendre polynomials at x that are needed
# by calling your code from prelab 
  p = eval_legendre(n,x)
  # initialize the sum to 0 
  pval = 0.0    
  for j in range(0,n+1):
      # make a function handle for evaluating phi_j(x)
      phi_j = lambda x: eval_last_legendre(j,x)
      # make a function handle for evaluating phi_j^2(x)*w(x)
      phi_j_sq = lambda x: ((phi_j(x))**2)*w(x)
      # use the quad function from scipy to evaluate normalizations
      norm_fac,err = quad(phi_j_sq,a,b)
      # make a function handle for phi_j(x)*f(x)*w(x)/norm_fac
      func_j = lambda x: phi_j(x)*f(x)*w(x)/norm_fac
      # use the quad function from scipy to evaluate coeffs
      aj,err = quad(func_j,a,b)
      # accumulate into pval
      pval = pval+aj*p[j] 
       
  return pval
#-------Lab Question Functions-------#
def PreLab():   
    N = 4
    x = 2

    p_test = eval_legendre(N,x)
    print(p_test)

def Problem3():
    #  function you want to approximate
    f = lambda x: math.exp(x)

# Interval of interest    
    a = -1
    b = 1
# weight function    
    w = lambda x: 1.

# order of approximation
    n = 2

#  Number of points you want to sample in [a,b]
    N = 1000
    xeval = np.linspace(a,b,N+1)
    pval = np.zeros(N+1)

    for kk in range(N+1):
      pval[kk] = eval_legendre_expansion(f,a,b,w,n,xeval[kk])
      
    ''' create vector with exact values'''
    fex = np.zeros(N+1)
    for kk in range(N+1):
        fex[kk] = f(xeval[kk])
        
    plt.figure()    
    plt.plot(xeval,fex,'ro-', label= 'f(x)')
    plt.plot(xeval,pval,'bs--',label= 'Expansion') 
    plt.legend()
    plt.title("Approximation for e^x")
    plt.show()    
    
    err = abs(pval-fex)
    plt.semilogy(xeval,err,'ro--',label='error')
    plt.title("Error for e^x")
    plt.legend()
    plt.show()

    #  function you want to approximate
    f = lambda x: 1/(1+x**2)

# Interval of interest    
    a = -1
    b = 1
# weight function    
    w = lambda x: 1.

# order of approximation
    n = 2

#  Number of points you want to sample in [a,b]
    N = 1000
    xeval = np.linspace(a,b,N+1)
    pval = np.zeros(N+1)

    for kk in range(N+1):
      pval[kk] = eval_legendre_expansion(f,a,b,w,n,xeval[kk])
      
    ''' create vector with exact values'''
    fex = np.zeros(N+1)
    for kk in range(N+1):
        fex[kk] = f(xeval[kk])
        
    plt.figure()    
    plt.plot(xeval,fex,'ro-', label= 'f(x)')
    plt.plot(xeval,pval,'bs--',label= 'Expansion') 
    plt.legend()
    plt.title("Approximation for 1/(1+x^2)")
    plt.show()    
    
    err = abs(pval-fex)
    plt.semilogy(xeval,err,'ro--',label='error')
    plt.title("Error for 1/(1+x^2)")
    plt.legend()
    plt.show()



#-------Calling Question Functions-------#
#PreLab()

Problem3()