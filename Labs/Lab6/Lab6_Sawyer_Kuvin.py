##APPM 4600 - Lab 5 (In Lab Code Attempt)
#Author: Sawyer Kuvin
#Date: 2/18/25
#Questions attempted/completed: 2.1, 2.2, 3.1, 3.2, ~3.3
# (~ partially completed/attempted) 

import numpy as np
from numpy import linalg as lg
import matplotlib.pyplot as plt
import scipy.special as sp
from scipy.linalg import lu_factor, lu_solve

disp = False

#-------Calculation functions-------#
def Forward_Difference(f,s,h):
    fprime = (f(s+h)-f(s))/h
    return(fprime)

def Centered_Difference(f,s,h):
    fprime = (f(s+h)-f(s-h))/(2*h)
    return(fprime)

def newton_method_nd(f,Jf,x0,tol,nmax,verb=False):

    # Initialize arrays and function value
    xn = x0; #initial guess
    rn = x0; #list of iterates
    Fn = f(xn); #function value vector
    n=0;
    nf=1; nJ=0; #function and Jacobian evals
    npn=1;

    if disp:
        print("|--n--|----xn----|---|f(xn)|---|");

    while npn>tol and n<=nmax:
        # compute n x n Jacobian matrix
        Jn = Jf(xn);
        nJ+=1;

        if disp:
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

    if disp:
        if np.linalg.norm(Fn)>tol:
            print("Newton method failed to converge, n=%d, |F(xn)|=%1.1e\n" % (nmax,np.linalg.norm(Fn)));
        else:
            print("Newton method converged, n=%d, |F(xn)|=%1.1e\n" % (n,np.linalg.norm(Fn)));

    return (r,rn,nf,nJ);

def slacker_newton_method_nd(f,Jf,x0,tol,nmax,verb=False):

    # Initialize arrays and function value
    xn = x0; #initial guess
    rn = x0; #list of iterates
    Fn = f(xn); #function value vector
    n=0;
    nf=1; nJ=0; #function and Jacobian evals
    npn=1;

    if (len(x0)<100):
        if (np.linalg.cond(Jf(x0)) > 1e16):
            print("Error: matrix too close to singular");
            print("Newton method failed to converge, n=%d, |F(xn)|=%1.1e\n" % (nmax,np.linalg.norm(Fn)));
            r=x0;
            return (r,rn,nf,nJ);

    if disp:
        print("|--n--|----xn----|---|f(xn)|---|");

    count = 0
    while npn>tol and n<=nmax:
        # compute n x n Jacobian matrix
        Jn = Jf(x0)
        if count == 5:
            Jn = Jf(xn);
            nJ+=1;
            count = 0

        if disp:
            print("|--%d--|%1.7f|%1.15f|" %(n,np.linalg.norm(xn),np.linalg.norm(Fn)));

        # Newton step (we could check whether Jn is close to singular here)
        pn = -np.linalg.solve(Jn,Fn);
        xn = xn + pn;
        npn = np.linalg.norm(pn); #size of Newton step

        n+=1;
        rn = np.vstack((rn,xn));
        Fn = f(xn);
        nf+=1;
        count += 1

    r=xn;

    if disp:
        if np.linalg.norm(Fn)>tol:
            print("Slacker Newton method failed to converge, n=%d, |F(xn)|=%1.1e\n" % (nmax,np.linalg.norm(Fn)));
        else:
            print("Slacker Newton method converged, n=%d, |F(xn)|=%1.1e\n" % (n,np.linalg.norm(Fn)));

    return (r,rn,nf,nJ);


def lazy_newton_method_nd(f,Jf,x0,tol,nmax,verb=False):

    # Initialize arrays and function value
    xn = x0; #initial guess
    rn = x0; #list of iterates
    Fn = f(xn); #function value vector
    # compute n x n Jacobian matrix (ONLY ONCE)
    Jn = Jf(xn);

    # Use pivoted LU factorization to solve systems for Jf. Makes lusolve O(n^2)
    lu, piv = lu_factor(Jn);

    n=0;
    nf=1; nJ=1; #function and Jacobian evals
    npn=1;

    if disp:
        print("|--n--|----xn----|---|f(xn)|---|");

    while npn>tol and n<=nmax:

        if disp:
            print("|--%d--|%1.7f|%1.12f|" %(n,np.linalg.norm(xn),np.linalg.norm(Fn)));

        # Newton step (we could check whether Jn is close to singular here)
        pn = -lu_solve((lu, piv), Fn); #We use lu solve instead of pn = -np.linalg.solve(Jn,Fn);
        xn = xn + pn;
        npn = np.linalg.norm(pn); #size of Newton step

        n+=1;
        rn = np.vstack((rn,xn));
        Fn = f(xn);
        nf+=1;

    r=xn;

    if disp:
        if np.linalg.norm(Fn)>tol:
            print("Lazy Newton method failed to converge, n=%d, |F(xn)|=%1.1e\n" % (nmax,np.linalg.norm(Fn)));
        else:
            print("Lazy Newton method converged, n=%d, |F(xn)|=%1.1e\n" % (n,np.linalg.norm(Fn)));

    return (r,rn,nf,nJ);

#-------Questions for Lab-------#
def BeforeLab():
    h = 0.01*2.0**(-np.arange(0,10))
    
    if disp == True:
        print(h)

    f = lambda x: np.cos(x)
    x0 = np.pi/2

    fprimeForward = Forward_Difference(f,x0,h)
   # print(fprimeForward)

    errorForward = fprimeForward-x0
    ForwardX = []
    ForwardY = []
    for j in range(len(errorForward)-1):
        ForwardX = np.append(ForwardX,errorForward[j])
        ForwardY = np.append(ForwardY,errorForward[j+1])
    alphaForward=alphaCenter = (ForwardY[3]-ForwardY[1])/(ForwardX[3]-ForwardX[1])
    plt.figure()
    plt.plot(ForwardX,ForwardY)
    plt.title("Forward Order Estimation")
    plt.show()
    print("The order of the Forward Difference is ", alphaForward)

    fprimeCentered = Centered_Difference(f,x0,h)
  #  print(fprimeCentered)

    errorCentered = fprimeCentered-x0
    CenterX = []
    CenterY = []
    for j in range(len(errorCentered)-1):
        CenterX = np.append(CenterX,errorCentered[j])
        CenterY = np.append(CenterY,errorCentered[j+1])
    alphaCenter = (CenterY[3]-CenterY[1])/(CenterX[3]-CenterX[1])
    plt.figure()
    plt.plot(CenterX,CenterY)
    plt.title("Centered Order Estimation")
    plt.show()
    print("The order of the Centered Difference is ", alphaCenter)

def Problem3():
    #Using original function from given code
    def F(x):
        return np.array([(x[0]-4)**2+2*(x[1]-2)**2-32 , x[1]*(x[0]-2)-16 ]);
    def JF(x):
        return np.array([[2*(x[0]-4),4*(x[1]-2)],[x[1],(x[0]-2)]]);

    x0 = np.array([4.0,4.0]); tol=1e-14; nmax=100;
    (rN,rnN,nfN,nJN) = newton_method_nd(F,JF,x0,tol,nmax,True);
    print("It took the Newton Method: ", nfN, " iterations and ", nJN, " Jacobians")

    (rSN,rnSN,nfSN,nJSN) = slacker_newton_method_nd(F,JF,x0,tol,nmax,True);
    print("It took the Slacker Method: ", nfSN, " iterations and ", nJSN, " Jacobians")

    (rLN,rnLN,nfLN,nJLN) = lazy_newton_method_nd(F,JF,x0,tol,nmax,True);
    print("It took the Lazy Method: ", nfLN, " iterations and ", nJLN, " Jacobians")

    def F2(x):
        return np.array([4*(x[0]**2) + x[1]**2-4 , x[0]+x[1]-np.sin(x[0]-x[1])]);
    def JF2(x):
        return np.array([[8*x[0],4*x[1]],[1-np.cos(x[0]-x[1]),1+np.cos(x[0]-x[1])]]);

    x02 = np.array([1,0])
    tol = 1e-10
    (rLN,rnLN,nfLN,nJLN) = lazy_newton_method_nd(F2,JF2,x02,tol,nmax);
    print("Lazy Method converged to: ", rLN, "in ", nfLN, " iterations")
    nmax = 100




#-------Running Lab Questions-------#
#BeforeLab()

Problem3()