##APPM 4600 - Homework 5
#Name: Sawyer Kuvin
#Date Started: 2/18/25
#Last Edit: 2/21/25

import numpy as np
from numpy import linalg as lg
import matplotlib.pyplot as plt
import scipy.special as sp
from scipy.linalg import lu_factor, lu_solve

#-------Calculation Functions-------#
def mod_lazy_newton_method_nd(f,x0,tol,nmax,verb=False):

    # Initialize arrays and function value
    xn = x0; #initial guess
    rn = x0; #list of iterates
    Fn = f(xn); #function value vector
    # compute n x n Jacobian matrix (ONLY ONCE)
    Jn_inv = np.array([[1/6,1/18],[0,1/6]])

    n=0;
    nf=1; nJ=1; #function and Jacobian evals
    npn=1;

    if verb:
        print("|--n--|----xn----|---|f(xn)|---|");

    while npn>tol and n<=nmax:

        if verb:
            print("|--%d--|%1.7f|%1.12f|" %(n,np.linalg.norm(xn),np.linalg.norm(Fn)));

        # Newton step (we could check whether Jn is close to singular here)
        pn = -Jn_inv @ Fn; #We use lu solve instead of pn = -np.linalg.solve(Jn,Fn);
        xn = xn + pn;
        npn = np.linalg.norm(pn); #size of Newton step

        n+=1;
        rn = np.vstack((rn,xn));
        Fn = f(xn);
        nf+=1;

    r=xn;

    if verb:
        if np.linalg.norm(Fn)>tol:
            print("Lazy Newton method failed to converge, n=%d, |F(xn)|=%1.1e\n" % (nmax,np.linalg.norm(Fn)));
        else:
            print("Lazy Newton method converged, n=%d, |F(xn)|=%1.1e\n" % (n,np.linalg.norm(Fn)));

    return (r,rn,nf,nJ);


def newton_method_nd(f,Jf,x0,tol,nmax,verb=False):

    # Initialize arrays and function value
    xn = x0; #initial guess
    rn = x0; #list of iterates
    Fn = f(xn); #function value vector
    n=0;
    nf=1; nJ=0; #function and Jacobian evals
    npn=1;

    if verb:
        print("|--n--|----xn----|---|f(xn)|---|");

    while npn>tol and n<=nmax:
        # compute n x n Jacobian matrix
        Jn = Jf(xn);
        nJ+=1;

        if verb:
            print("|--%d--|%1.7f|%1.12f|" %(n,np.linalg.norm(xn),np.linalg.norm(Fn)));

        # Newton step (we could check whether Jn is close to singular here)
        pn = -np.linalg.solve(Jn,Fn);
        xn = xn + pn;
        npn = np.linalg.norm(pn); #size of Newton step

        n+=1;
        rn = np.vstack((rn,xn));
        Fn = f(xn);
        nf+=1;

    r=xn;

    if verb:
        if np.linalg.norm(Fn)>tol:
            print("Newton method failed to converge, n=%d, |F(xn)|=%1.1e\n" % (nmax,np.linalg.norm(Fn)));
        else:
            print("Newton method converged, n=%d, |F(xn)|=%1.1e\n" % (n,np.linalg.norm(Fn)));

    return (r,rn,nf,nJ);

def plane_intersect_solve(f,fprime,x0,tol,nmax):
     # Initialize arrays and function value
    xn = x0; #initial guess
    rn = x0; #list of iterates
    Fn = f(xn); #function value vector
    n=0;
    nf=1; nP=0; #function and prime evals
    npn=1;
    Nmax = nmax

    while npn>tol and n<Nmax:
        fprimen = fprime(xn)
        Fn = f(xn)
        nP+=1
        nf+=1

        dx = Fn / (fprimen[0]**2+fprimen[1]**2+fprimen[2]**2)
        dy = Fn / (fprimen[0]**2+fprimen[1]**2+fprimen[2]**2)
        dz = Fn / (fprimen[0]**2+fprimen[1]**2+fprimen[2]**2)
        xn = np.array([(xn[0]-dx*fprimen[0]),(xn[1]-dy*fprimen[1]),(xn[2]-dz*fprimen[2])])
        rn = np.vstack((rn,xn))
        n+=1
        #print(n)
    
    r = xn
    return (r,rn,nf,nf)




#-------Problem Functions-------#
def Problem1():
    #Part a)
    def F(x):
        return np.array([3*(x[0]**2)-x[1]**2, 3*(x[0]*(x[1]**2))-(x[0]**3)-1])
    x0 = np.array([1,1])
    tol = 1e-10
    nmax = 100
    (rI,rnI,nfI,_) = mod_lazy_newton_method_nd(F,x0,tol,nmax)
    print("The given method converges to: ", rI, " and it took ", nfI, " iterations.")
    
    root = rI
    initialError = abs(rnI - root)

    initialX = []
    initialY = []
    for j in range(len(rnI)-2):
        initialX = np.append(initialX, initialError[j])
        initialY = np.append(initialY, initialError[j+1])
    initialOrder = (np.log10(initialY[-2])-np.log10(initialY[2]))/(np.log10(initialX[-2])-np.log10(initialX[2]))
    plt.figure()
    plt.loglog(initialX,initialY)
    plt.grid(True)
    plt.xlabel("k Iteration Error")
    plt.ylabel("k+1 Iteration Error")
    plt.title("1.a - Convergance of Intial Matrix Method")
    plt.show()

    print("The order of the initial method was: ", initialOrder)

    #Part b)
    def JF(x):
        return np.array([[6*x[0],-2*x[1]],[3*(x[1]**2)-3*(x[0]**2),6*x[0]*x[1]]])
    
    (rN,rnN,nfN,_) = newton_method_nd(F,JF,x0,tol,nmax,verb=False)
    print("The Newton method converges to: ", rN, " and it took ", nfN, " iterations.")
    
    root = rN
    NewtonError = abs(rnN - root)

    NewtonX = []
    NewtonY = []
    for j in range(len(rnN)-2):
        NewtonX = np.append(NewtonX, NewtonError[j])
        NewtonY = np.append(NewtonY, NewtonError[j+1])
    NewtonOrder = (np.log10(NewtonY[-3])-np.log10(NewtonY[3]))/(np.log10(NewtonX[-3])-np.log10(NewtonX[3]))
    plt.figure()
    plt.loglog(initialX,initialY)
    plt.loglog(NewtonX[0:-1],NewtonY[0:-1])
    plt.grid(True)
    plt.xlabel("k Iteration Error")
    plt.ylabel("k+1 Iteration Error")
    plt.legend(["Initial Method (Order ~ 1)", "Newton Method (Order ~ 2)"])
    plt.title("1.c - Convergance of Two Methods")
    plt.show()

    print("The order of the Newton method was: ", NewtonOrder)

def Problem2():
    x_vec = np.linspace(-2,2,20)
    f_left = lambda x: -1-x
    f_right = lambda x: 1-x
    g_left = lambda x: -1+x
    g_right = lambda x: 1+x

    plt.figure()
    plt.plot(x_vec, f_left(x_vec), color="blue")
    plt.plot(x_vec, g_left(x_vec), color="red")
    plt.plot(x_vec, f_right(x_vec), color="blue",label='_nolegend_')
    plt.plot(x_vec, g_right(x_vec), color="red",label='_nolegend_')
    plt.fill_between(x_vec,f_left(x_vec),f_right(x_vec), color="blue", alpha=0.5, label='_nolegend_')
    plt.fill_between(x_vec,g_left(x_vec),g_right(x_vec), color="red", alpha=0.5, label='_nolegend_')
    plt.legend(["Bounds of 1st Equation", "Bounds of 2nd Equation"])
    plt.ylim([-2,2])
    plt.xlabel("X Points")
    plt.ylabel("Y Points")
    plt.title("Domain w/ Guaranteed Fixed Point Convergence")
    plt.show()


def Problem3():
    def F(x):
        return x[0]**2+4*(x[1]**2)+4*(x[2]**2)-16
    def Fprime(x):
        return np.array([2*x[0],8*x[1],8*x[2]])
    nmax = 5
    tol = 1e-10

    x0 = np.array([1,1,1])

    (root,rN,_,_) = plane_intersect_solve(F,Fprime,x0,tol,nmax)

    print("The given method finds an intersection at: ", root)
    Error = abs(rN-root)

    ErrorX = []
    ErrorY = []
    for j in range(len(rN)-2):
        ErrorX = np.append(ErrorX, Error[j])
        ErrorY = np.append(ErrorY, Error[j+1])
    ErrorOrder = (np.log10(ErrorY[-1])-np.log10(ErrorY[1]))/(np.log10(ErrorX[-1])-np.log10(ErrorX[1]))
    plt.figure()
    plt.loglog(ErrorX,ErrorY)
    plt.loglog(ErrorX[0:-1],ErrorY[0:-1])
    plt.grid(True)
    plt.xlabel("k Iteration Error")
    plt.ylabel("k+1 Iteration Error")
    plt.title("3.b - Convergance of Intersection Iterations")
    plt.show()

    print("The order of the method was: ", ErrorOrder)
    

#-------Calling Problem Functions-------#
#Problem 1 call
Problem1()

#Problem 2 call
Problem2()

#Problem 3 call
Problem3()



