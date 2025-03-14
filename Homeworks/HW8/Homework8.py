##APPM 4600 - Homework 8
#Author: Sawyer Kuvin
#Date: 3/14/25


#Importing Commonly Used Libraries
import numpy as np
from numpy import linalg as lg
import matplotlib.pyplot as plt
import scipy.special as sp
from scipy.linalg import lu_factor, lu_solve, inv

#--------Interpretation Method Functions-------#
def eval_lagrange(xeval,xint,yint,N): #Given by Professor

    lj = np.ones(N+1)
    for count in range(N+1):
       for jj in range(N+1):
           if (jj != count):
              lj[count] = lj[count]*(xeval - xint[jj])/(xint[count]-xint[jj])

    yeval = 0.
    
    for jj in range(N+1):
       yeval = yeval + yint[jj]*lj[jj]
  
    return(yeval)

def eval_hermite(xeval,xint,yint,yprimeint,N):
    lj = np.ones(N+1)
    for count in range(N+1):
       for jj in range(N+1):
           if (jj != count):
              lj[count] = lj[count]*(xeval - xint[jj])/(xint[count]-xint[jj])

    lprimej = np.zeros(N+1)

    for count in range(N+1):
        for jj in range(N+1):
            if (jj != count):
                lprimej[count] = lprimej[count]+ 1./(xint[count]-xint[jj])
    
    yeval = 0.
    for jj in range(N+1):
        Hk = (1.-(2.*(xeval-xint[jj])*lprimej[jj]))*lj[jj]**2
        Kk = (xeval-xint[jj])*lj[jj]**2

        yeval = yeval + yint[jj]*Hk+ yprimeint[jj]*Kk

    return(yeval)

def create_natural_spline(yint,xint,N): #Given by Professor
    # create the right hand side for the linear system
    b = np.zeros(N+1)
    # vector values
    h = np.zeros(N+1)
    for i in range(1,N):
        hi = xint[i]-xint[i-1]
        hip = xint[i+1] - xint[i]
        b[i] = (yint[i+1]-yint[i])/hip - (yint[i]-yint[i-1])/hi
        h[i-1] = hi
        h[i] = hip
    # create matrix so you can solve for the M values
    # This is made by filling one row at a time
    A = np.zeros((N+1,N+1))
    A[0][0] = 1.0
    for j in range(1,N):
        A[j][j-1] = h[j-1]/6
        A[j][j] = (h[j]+h[j-1])/3
        A[j][j+1] = h[j]/6
    A[N][N] = 1
    Ainv = inv(A)
    M = Ainv.dot(b)
    # Create the linear coefficients
    C = np.zeros(N)
    D = np.zeros(N)
    for j in range(N):
        C[j] = yint[j]/h[j]-h[j]*M[j]/6
        D[j] = yint[j+1]/h[j]-h[j]*M[j+1]/6
    return(M,C,D)

def create_clamped_spline(yint,xint,f_prime_a,f_prime_b,N):
    b = np.zeros(N+1)
    h = np.zeros(N+1)
    for i in range(1,N):
        hi = xint[i]-xint[i-1]
        hip = xint[i+1] - xint[i]
        b[i] = (yint[i+1]-yint[i])/hip - (yint[i]-yint[i-1])/hi
        h[i-1] = hi
        h[i] = hip
    b[0] = (yint[1]-yint[0])/h[0] - f_prime_a
    b[N] = f_prime_b - (yint[N]-yint[N-1])/h[N-1]

    A = np.zeros((N+1,N+1))
    A[0][0] = 2*h[0]
    A[1][0] = h[0]
    A[0][1] = h[0]
    for j in range(1,N):
        A[j][j-1] = h[j-1]/6
        A[j][j] = (h[j]+h[j-1])/3
        A[j][j+1] = h[j]/6
    A[N][N] = 2*h[N-1]
    A[N-1][N] = h[N-1]
    A[N][N-1] = h[N-1]

    Ainv = inv(A)
    M = Ainv.dot(b)
    # Create the linear coefficients
    C = np.zeros(N)
    D = np.zeros(N)
    for j in range(N):
        C[j] = yint[j]/h[j]-h[j]*M[j]/6
        D[j] = yint[j+1]/h[j]-h[j]*M[j+1]/6
    return(M,C,D)

def create_periodic_spline(yint,xint,N):
    # create the right hand side for the linear system
    b = np.zeros(N+1)
    # vector values
    h = np.zeros(N+1)
    for i in range(1,N):
        hi = xint[i]-xint[i-1]
        hip = xint[i+1] - xint[i]
        b[i] = (yint[i+1]-yint[i])/hip - (yint[i]-yint[i-1])/hi
        h[i-1] = hi
        h[i] = hip

    A = np.zeros((N+1,N+1))
    A[0][0] = 1.0
    A[0][N] = h[0]/6
    for j in range(1,N):
        A[j][j-1] = h[j-1]/6
        A[j][j] = (h[j]+h[j-1])/3
        A[j][j+1] = h[j]/6
    A[N][N] = 1.0
    A[N][0] = h[N]/6
    Ainv = inv(A)
    M = Ainv.dot(b)
    # Create the linear coefficients
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
    yeval = (Mi*(xip-xeval)**3 +(xeval-xi)**3*Mip)/(6*hi) \
            + C*(xip-xeval) + D*(xeval-xi)
    return yeval

def eval_cubic_spline(xeval,Neval,xint,Nint,M,C,D):
    yeval = np.zeros(Neval+1)
    for j in range(Nint):
        '''find indices of xeval in interval (xint(jint),xint(jint+1))'''
        '''let ind denote the indices in the intervals'''
        atmp = xint[j]
        btmp= xint[j+1]
        # find indices of values of xeval in the interval
        ind= np.where((xeval >= atmp) & (xeval <= btmp))
        xloc = xeval[ind]
        # evaluate the spline
        yloc = eval_local_spline(xloc,atmp,btmp,M[j],M[j+1],C[j],D[j])
        # print('yloc = ', yloc)
        # copy into yeval
        yeval[ind] = yloc
    return(yeval)

#-------Problem Functions-------#
def Problem1():
    Nvec = [5,10,15,20]
    #Nvec = [5]
    print(Nvec)
    a = -5
    b = 5
    for i in range(len(Nvec)):
        N = Nvec[i]
        #Equispaced Nodes
        x_i = np.linspace(a,b,N+1)

        f_x = lambda x: 1 / (1+x**2)
        f_prime_x = lambda x: (-2*x)/((1+x**2)**2)
        y_i = f_x(x_i)
        y_prime_i = f_prime_x(x_i)

        Neval = 101
        x_test = np.linspace(a,b,Neval)
        y_real = f_x(x_test)

        #Lagrange
        yeval_l= np.zeros(Neval)
        y = np.zeros( (N+1, N+1) )
        for j in range(N):
            y[j][0]  = y_i[j]
        for kk in range(Neval):
            yeval_l[kk] = eval_lagrange(x_test[kk],x_i,y_i,N)

        #Hermite
        yeval_h= np.zeros(Neval)
        y = np.zeros( (N+1, N+1) )
        for j in range(N):
            y[j][0]  = y_i[j]
        for kk in range(Neval):
            yeval_h[kk] = eval_hermite(x_test[kk],x_i,y_i,y_prime_i,N)

        #Natural Cubic Spline
        (M_nat,C_nat,D_nat) = create_natural_spline(y_i,x_i,N)
        yeval_nat = eval_cubic_spline(x_test,(Neval-1),x_i,N,M_nat,C_nat,D_nat)


        #Clamped Cubic Spline
        f_prime_a = f_prime_x(a)
        f_prime_b = f_prime_x(b)
        (M_clam,C_clam,D_clam) = create_clamped_spline(y_i,x_i,f_prime_a,f_prime_b,N)
        yeval_clam = eval_cubic_spline(x_test,(Neval-1),x_i,N,M_clam,C_clam,D_clam)

        plt.figure()
        plt.plot(x_test,y_real,'k-', label='Original Function')
        plt.plot(x_test, yeval_l, 'r--', label='Lagrange Interpolation')
        plt.plot(x_test, yeval_h, 'b--', label='Hermite Interpolation')
        plt.plot(x_test, yeval_nat, 'g--', label='Natural Cubic Spline')
        plt.plot(x_test, yeval_clam, 'm--', label='Clamped Cubic Spline')
        plt.title("Different Interpolated Plots for N = "+ str(N))
        plt.grid(True)
        plt.grid('minor')
        plt.ylim([-0.25,1.25])
        plt.legend()
        plt.show()

        plt.figure()
        #plt.scatter(x_i,y_i,s=25, c='r')
        #err_v = abs(y_real - y_mon)python
        err_l = abs(y_real - yeval_l)
        err_h = abs(y_real - yeval_h)
        err_nat = abs(y_real - yeval_nat)
        err_clam = abs(y_real - yeval_clam)
        plt.semilogy(x_test,err_l,'ro--',label='Lagrange')
        plt.semilogy(x_test,err_h,'bo--',label='Hermite')
        plt.semilogy(x_test,err_nat,'go--',label='Natural Spline')
        plt.semilogy(x_test,err_clam,'mo--',label='Clamped Spline')
        plt.grid(True)
        plt.grid('minor')
        plt.title("Equispaced - Interpolation Error for N = " + str(N))
        plt.ylabel("Error between Model and Function")
        plt.xlabel("X points")
        plt.legend()
        plt.show()

def Problem2():
    Nvec = [5,10,15,20]
    print(Nvec)
    #Nvec = [5]
    a = -5
    b = 5
    for i in range(len(Nvec)):
        N = Nvec[i]
        #Chebushev Nodes
        x_i = np.zeros(N+2)
        for i in range(1,N+2):
            x_i[i] = np.cos(((2*i-1)*np.pi)/(2*N+2))
        
        x_i = np.flip(5*x_i)
        x_i = x_i[0:N+1]
        #print(x_i)

        f_x = lambda x: 1 / (1+x**2)
        f_prime_x = lambda x: (-2*x)/((1+x**2)**2)
        y_i = f_x(x_i)
        #print(y_i)
        y_prime_i = f_prime_x(x_i)

        Neval = 101
        x_test = np.linspace(a,b,Neval)
        y_real = f_x(x_test)

        #Lagrange
        yeval_l= np.zeros(Neval)
        y = np.zeros( (N+1, N+1) )
        for j in range(N):
            y[j][0]  = y_i[j]
        for kk in range(Neval):
            yeval_l[kk] = eval_lagrange(x_test[kk],x_i,y_i,N)

        #Hermite
        yeval_h= np.zeros(Neval)
        y = np.zeros( (N+1, N+1) )
        for j in range(N):
            y[j][0]  = y_i[j]
        for kk in range(Neval):
            yeval_h[kk] = eval_hermite(x_test[kk],x_i,y_i,y_prime_i,N)

        #Natural Cubic Spline
        (M_nat,C_nat,D_nat) = create_natural_spline(y_i,x_i,N)
        yeval_nat = eval_cubic_spline(x_test,(Neval-1),x_i,N,M_nat,C_nat,D_nat)


        #Clamped Cubic Spline
        f_prime_a = f_prime_x(a)
        f_prime_b = f_prime_x(b)
        (M_clam,C_clam,D_clam) = create_clamped_spline(y_i,x_i,f_prime_a,f_prime_b,N)
        yeval_clam = eval_cubic_spline(x_test,(Neval-1),x_i,N,M_clam,C_clam,D_clam)

        plt.figure()
        plt.plot(x_test,y_real,'k-', label='Original Function')
        plt.plot(x_test, yeval_l, 'r--', label='Lagrange Interpolation')
        plt.plot(x_test, yeval_h, 'b--', label='Hermite Interpolation')
        plt.plot(x_test, yeval_nat, 'g--', label='Natural Cubic Spline')
        plt.plot(x_test, yeval_clam, 'm--', label='Clamped Cubic Spline')
        plt.title("Chebyshev - Different Interpolated Plots for N = "+ str(N))
        plt.grid(True)
        plt.grid('minor')
        plt.ylim([-0.25,1.25])
        plt.legend()
        plt.show()

        plt.figure()
        #plt.scatter(x_i,y_i,s=25, c='r')
        #err_v = abs(y_real - y_mon)python
        err_l = abs(y_real - yeval_l)
        err_h = abs(y_real - yeval_h)
        err_nat = abs(y_real - yeval_nat)
        err_clam = abs(y_real - yeval_clam)
        plt.semilogy(x_test,err_l,'ro--',label='Lagrange')
        plt.semilogy(x_test,err_h,'bo--',label='Hermite')
        plt.semilogy(x_test,err_nat,'go--',label='Natural Spline')
        plt.semilogy(x_test,err_clam,'mo--',label='Clamped Spline')
        plt.grid(True)
        plt.grid('minor')
        plt.title("Chebyshev - Interpolation Error for N = " + str(N))
        plt.ylabel("Error between Model and Function")
        plt.xlabel("X points")
        plt.legend()
        plt.show()

def Problem3():
    N_vec = [5,10,13,27,40]
    for n in range(len(N_vec)):
        N = N_vec[n]
        a = 0
        b = 2 * np.pi

        x_i = np.linspace(a,b,N+1)

        f_x = lambda x: np.sin(10*x)
        y_i = f_x(x_i)

        Neval = 101
        x_test = np.linspace(a,b,Neval)
        y_real = f_x(x_test)

        (M_per,C_per,D_per) = create_periodic_spline(y_i,x_i,N)
        yeval_per = eval_cubic_spline(x_test,(Neval-1),x_i,N,M_per,C_per,D_per)

        plt.figure()
        plt.plot(x_test,y_real,'k-', label='Original Function')
        plt.plot(x_test, yeval_per, 'r--', label='Periodic Spline')
        plt.title("Periodic Interpolation Test for N = " + str(N))
        plt.grid(True)
        plt.grid('minor')
        plt.legend()
        plt.show()

        plt.figure()
        #plt.scatter(x_i,y_i,s=25, c='r')
        #err_v = abs(y_real - y_mon)python
        err_per = abs(y_real - yeval_per)
        plt.semilogy(x_test,err_per,'ro--',label='Periodic Spline')
        if N == 40:
            plt.title("Plot with No Aliasing: N = " + str(N))
        else:
            plt.title("Plot with Aliasing: N = " + str(N))
        plt.grid(True)
        plt.grid('minor')
        plt.show()
    


#-------Calling Problem Functions-------#
Problem1()

Problem2()

Problem3()