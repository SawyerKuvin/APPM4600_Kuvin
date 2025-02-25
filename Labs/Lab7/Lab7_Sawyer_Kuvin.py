##APPM 4600 - Lab 7 (In Lab Code Attempt)
#Author: Sawyer Kuvin
#Date: 2/25/25
#Questions attempted/completed: 2, ~3.1
# (~ partially completed/attempted) 

#Importing Commonly Used Libraries
import numpy as np
from numpy import linalg as lg
import matplotlib.pyplot as plt
import scipy.special as sp
from scipy.linalg import lu_factor, lu_solve


#--------Reused/Calculation Functions-------#
def  eval_lin_spline(xeval,Neval,a,b,f,Nint): #Given Function

    '''create the intervals for piecewise approximations'''
    xint = np.linspace(a,b,Nint+1)
   
    '''create vector to store the evaluation of the linear splines'''
    yeval = np.zeros(Neval) 
    
    for jint in range(Nint):
        '''find indices of xeval in interval (xint(jint),xint(jint+1))'''
        '''let ind denote the indices in the intervals'''
        '''let n denote the length of ind'''
        
        '''temporarily store your info for creating a line in the interval of 
         interest'''
        a1= xint(jint)
        fa1 = f(a1)
        b1 = xint(jint+1)
        fb1 = f(b1)
        
        for kk in range(n):
           '''use your line evaluator to evaluate the lines at each of the points 
           in the interval'''
           '''yeval(ind(kk)) = call your line evaluator at xeval(ind(kk)) with 
           the points (a1,fa1) and (b1,fb1)'''

def eval_lagrange(xeval,xint,yint,N):

    lj = np.ones(N)
    
    for count in range(N):
       for jj in range(N):
           if (jj != count):
              lj[count] = lj[count]*(xeval - xint[jj])/(xint[count]-xint[jj])

    yeval = 0.
    
    for jj in range(N):
       yeval = yeval + yint[jj]*lj[jj]
  
    return(yeval)

def dividedDiffTable(x, y, n):
 
    for i in range(1, n):
        for j in range(n - i):
            y[j][i] = ((y[j][i - 1] - y[j + 1][i - 1]) /
                                     (x[j] - x[i + j]));
    return y;
    
def evalDDpoly(xval, xint,y,N):
    ''' evaluate the polynomial terms'''
    ptmp = np.zeros(N+1)
    
    ptmp[0] = 1.
    for j in range(N):
      ptmp[j+1] = ptmp[j]*(xval-xint[j])
     
    '''evaluate the divided difference polynomial'''
    yeval = 0.
    for j in range(N):
       yeval = yeval + y[0][j]*ptmp[j]  

    return yeval

def  eval_monomial(xeval,coef,N,Neval):

    yeval = coef[0]*np.ones(Neval)
    
#    print('yeval = ', yeval)
    
    for j in range(1,N):
      for i in range(Neval):
#        print('yeval[i] = ', yeval[i])
#        print('a[j] = ', a[j])
#        print('i = ', i)
#        print('xeval[i] = ', xeval[i])
        yeval[i] = yeval[i] + coef[j]*xeval[i]**j

    return yeval

def vandermonde(x_vec): #Coded myself
    #x_vec - Vector of x values from original function
    #N - length of y vector of values
    N = len(x_vec)
    V = np.zeros((N,N))
    for i in range(N):
        for j in range(N):
            V[i,j] = x_vec[i]**j

    return(V)

def interp_coeffs(V,y_vec): #Coded myself
    a = np.linalg.solve(V,y_vec)
    return(a)

def monomial_expansion(x_vec,y_vec,x_0=-7000): #Coded myself
    V = vandermonde(x_vec)
    a = interp_coeffs(V,y_vec)

    if x_0==-7000:
        return(V,a)
    else:
        p_0 = 0
        N = len(x_vec)
        for k in range(N):
            p_0 = p_0 + a[k]*(x_0**k)
        return(p_0,V,a)
    
#-------Problem Functions--------#
def Problem3_1():
    #Monomial Expansion
    f_x = lambda x: 1/(1+(10*x)**2)
    N = 18
    '''j = np.arange(0,N,1)
    h = 2/(N-1)
    x_j = -1 + ((j-1)/h)'''
    x_j = np.linspace(-1,1,N)
    print(x_j)
    y_j = f_x(x_j)
    [_,a_j] = monomial_expansion(x_j,y_j)
    print(a_j)

    Neval = 1000
    x_test = np.linspace(-1,1,Neval)
    y_mon = eval_monomial(x_test,a_j,N,Neval)
    y_real = f_x(x_test)

    plt.figure()
    plt.scatter(x_j,y_j)
    plt.plot(x_test,y_mon)
    plt.plot(x_test,y_real)
    plt.legend("Monomial Estimation", "Actual Function")
    plt.title("Monomial Interpolation for N = 18")
    plt.grid(True)
    plt.show()

    #Langrange Polynomials
    
    #Failed attempt at coding by myself
    '''xeval = np.linspace(-1,1,Neval)
    yeval_l= np.zeros(Neval)
    yeval_dd = np.zeros(Neval)

    y = np.zeros( (N, N) )

    for j in range(N):
       y[j][0]  = y_j[j]
    
    for kk in range(Neval):
       yeval_l[kk] = eval_lagrange(xeval[kk],x_j,y_j,N)
       yeval_dd[kk] = evalDDpoly(xeval[kk],x_j,y,N)

    fex = f_x(xeval)
       

    plt.figure()    
    plt.plot(xeval,fex,'ro-')
    plt.plot(xeval,yeval_l,'bs--') 
    plt.plot(xeval,yeval_dd,'c.--')
    plt.legend()
    plt.show()

    plt.figure() 
    err_l = abs(yeval_l-fex)
    err_dd = abs(yeval_dd-fex)
    plt.semilogy(xeval,err_l,'ro--',label='lagrange')
    plt.semilogy(xeval,err_dd,'bs--',label='Newton DD')
    plt.legend()
    plt.show()'''

    ''' create equispaced interpolation nodes'''
    xint = np.linspace(-1,1,N+1)
    
    ''' create interpolation data'''
    yint = f_x(xint)
    
    ''' create points for evaluating the Lagrange interpolating polynomial'''
    Neval = 1000
    xeval = np.linspace(-1,1,Neval+1)
    yeval_l= np.zeros(Neval+1)
    yeval_dd = np.zeros(Neval+1)
  
    '''Initialize and populate the first columns of the 
     divided difference matrix. We will pass the x vector'''
    y = np.zeros( (N+1, N+1) )
     
    for j in range(N+1):
       y[j][0]  = yint[j]

    y = dividedDiffTable(xint, y, N+1)
    ''' evaluate lagrange poly '''
    for kk in range(Neval+1):
       yeval_l[kk] = eval_lagrange(xeval[kk],xint,yint,N)
       yeval_dd[kk] = evalDDpoly(xeval[kk],xint,y,N)
          


    ''' create vector with exact values'''
    fex = f_x(xeval)
       

    plt.figure()    
    plt.plot(xeval,fex,'ro-')
    plt.plot(xeval,yeval_l,'bs--') 
    plt.plot(xeval,yeval_dd,'c.--')
    plt.legend()
    plt.title("Langrange and Newton Interpolation for N = 18")

    plt.figure() 
    err_l = abs(yeval_l-fex)
    err_dd = abs(yeval_dd-fex)
    plt.semilogy(xeval,err_l,'ro--',label='lagrange')
    plt.semilogy(xeval,err_dd,'bs--',label='Newton DD')
    plt.legend()
    plt.title("Error in Langrange and Newton Interpolation for N = 18")
    plt.show()

#-------Calling Problems-------#
Problem3_1()

#Note: Wrote code for Vandermonde by myself and coded up the langrange
# and newton dd with the code given by the proffessor.

