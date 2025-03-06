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

#-------Given/Made for Lab Functions-------#
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

def eval_lagrange(xeval,xint,yint,N): #Given by Professor

    lj = np.ones(N)
    
    for count in range(N):
       for jj in range(N):
           if (jj != count):
              lj[count] = lj[count]*(xeval - xint[jj])/(xint[count]-xint[jj])

    yeval = 0.
    
    for jj in range(N):
       yeval = yeval + yint[jj]*lj[jj]
  
    return(yeval)


#-------Problem Functions-------#
def Problem1():
    N_vec = [2,3,4,5,16,17,18,19]
    for i in range(len(N_vec)):
        N = N_vec[i]
        h = 2/(N-1)
        x_i = x_n = -1
        for i in range(N-1):
            x_n = -1 + (i+1)*h
            x_i = np.vstack((x_i,x_n))
        #print(x_i)

        f_x = lambda x: 1/(1+((10*x)**2))
        y_i = f_x(x_i)
        [_,c_n] = monomial_expansion(x_i,y_i)
        
        Neval = 1001
        x_test = np.linspace(-1,1,Neval)
        y_mon = eval_monomial(x_test,c_n,N,Neval)
        y_real = f_x(x_test)

        if abs(max(y_mon)) < 100:
            plt.figure()
            plt.scatter(x_i,y_i,s=25, c='r')
            plt.plot(x_test,y_mon)
            plt.plot(x_test,y_real)
            plt.grid('minor')
            plt.legend(["Interpolation Points","Vandermonde Estimation", "Actual Function"])
            plt.title("Vandermonde Interpolation for N = " + str(N))
            plt.grid(True)
            plt.show()
        elif abs(max(y_mon)) > 100:
            print("Exceeds Tolerance")
            plt.figure()
            plt.scatter(x_i,y_i,s=25, c='r')
            plt.plot(x_test,y_mon)
            plt.plot(x_test,y_real)
            plt.grid('minor')
            plt.legend(["Interpolation Points","Vandermonde Estimation", "Actual Function"])
            plt.title("Vandermonde Interpolation for N = " + str(N))
            plt.grid(True)
            plt.show() 

def Problem2():
    N_vec = [2,3,4,5,16,17,18,19,100,150,200]
    for i in range(len(N_vec)):
        N = N_vec[i]
        h = 2/(N-1)
        x_i = x_n = -1
        for i in range(N-1):
            x_n = -1 + (i+1)*h
            x_i = np.vstack((x_i,x_n))
        #print(x_i)

        f_x = lambda x: 1/(1+((10*x)**2))
        y_i = f_x(x_i)
        
        Neval = 1001
        x_test = np.linspace(-1,1,Neval)
        y_real = f_x(x_test)

    #Code given by professor for Lagrange Interpolation
        yeval_l= np.zeros(Neval)
        y = np.zeros( (N+1, N+1) )
        for j in range(N):
            y[j][0]  = y_i[j]
        for kk in range(Neval):
            yeval_l[kk] = eval_lagrange(x_test[kk],x_i,y_i,N)

        if N < 100:
            plt.figure()
            #plt.scatter(x_i,y_i,s=25, c='r')
            #err_v = abs(y_real - y_mon)
            err_l = abs(y_real - yeval_l)
            #plt.semilogy(x_test,err_v,'bs--',label='vandermonde')
            plt.semilogy(x_test,err_l,'ro--',label='lagrange')
            plt.grid('minor')
            plt.title("Lagrange Interpolation Error for N = " + str(N))
            plt.grid(True)
            plt.show()

        if N >= 100:
            plt.figure()
            #plt.scatter(x_i,y_i,s=25, c='r')
            #err_v = abs(y_real - y_mon)
            err_l = abs(y_real - yeval_l)
            #plt.semilogy(x_test,err_v,'bs--',label='vandermonde')
            plt.semilogy(x_test,err_l,'ro--',label='lagrange')
            plt.grid('minor')
            plt.title("Lagrange Interpolation Error for N = " + str(N))
            plt.grid(True)
            plt.show()
        

def Problem3():
    N_vec = [2,3,4,5,16,17,18,19]
    for i in range(len(N_vec)):
        N = N_vec[i]
        h = 2/(N-1)
        x_i = -1
        for i in range(N-1):
            x_n = -1 + (i+1)*h
            x_i = np.vstack((x_i,x_n))
        #print(x_i)

        x_j = np.cos(np.pi/(2*N))
        for j in range(N-1):
            x_n = np.cos((((2*(j+2))-1)*np.pi)/(2*N))
            x_j = np.vstack((x_j,x_n))

        #print(x_j)

        f_x = lambda x: 1/(1+((10*x)**2))
        y_i = f_x(x_i)
        y_j = f_x(x_j)
        
        Neval = 1001
        x_test = np.linspace(-1,1,Neval)
        y_real = f_x(x_test)

    #Code given by professor for Lagrange Interpolation
        yeval_l_e= np.zeros(Neval)
        y = np.zeros( (N+1, N+1) )
        for j in range(N):
            y[j]  = y_i[j]
        for kk in range(Neval):
            yeval_l_e[kk] = eval_lagrange(x_test[kk],x_i,y_i,N)

        yeval_l_c= np.zeros(Neval)
        y = np.zeros( (N+1, N+1) )
        for j in range(N):
            y[j]  = y_j[j]
        for kk in range(Neval):
            yeval_l_c[kk] = eval_lagrange(x_test[kk],x_j,y_j,N)

        '''plt.figure()
        #plt.scatter(x_i,y_i,s=25, c='r')
        #err_v = abs(y_real - y_mon)
        err_l_e = abs(y_real - yeval_l_e)
        err_l_c = abs(y_real - yeval_l_c)
        plt.semilogy(x_test,err_l_e,'ro--',label='Equispaced Nodes')
        plt.semilogy(x_test,err_l_c,'bs--',label='Chebushev Nodes')
        plt.grid(True)
        plt.grid('minor')
        plt.legend()
        plt.title("Chebushev Error Comparison for N = " + str(N))
        plt.show()'''

        if N == 16:
            plt.figure()
            plt.plot(x_test,y_real, label='Original Function')
            plt.plot(x_test, yeval_l_e, label='Equispaced Nodes')
            plt.plot(x_test, yeval_l_c, label='Chebushev Nodes')
            plt.grid(True)
            plt.grid('minor')
            plt.legend()
            plt.title("Different Node Type Interpolation")
            plt.show()

#-------Calling Problem Functions-------#
#Problem 1
Problem1()

#Problem 2
Problem2()

#Problem 3
Problem3()