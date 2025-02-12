##APPM 4600 - Homework 4
#Name: Sawyer Kuvin
#Date Started: 2/11/25
#Last Edit: 2/11/25

import numpy as np
import matplotlib.pyplot as plt
import scipy.special as sp

#-------Given Functions-------#
#Bisection Method
def bisection(f,a,b,tol):
    
#    Inputs:
#     f,a,b       - function and endpoints of initial interval
#      tol  - bisection stops when interval length < tol
#      num  - number of iterations

#    Returns:
#      astar - approximation of root
#      ier   - error message
#            - ier = 1 => Failed
#            - ier = 0 == success

#     first verify there is a root we can find in the interval 
    num = 0
    fa = f(a)
    fb = f(b);
    if (fa*fb>0):
       ier = 1
       astar = a
       return [astar, ier,num]

#   verify end points are not a root 
    if (fa == 0):
      astar = a
      ier =0
      return [astar, ier,num]

    if (fb ==0):
      astar = b
      ier = 0
      return [astar, ier,num]

    count = 0
    d = 0.5*(a+b)
    while (abs(d-a)> tol):
      num = num + 1
      fd = f(d)
      if (fd ==0):
        astar = d
        ier = 0
        return [astar, ier, num]
      if (fa*fd<0):
         b = d
      else: 
        a = d
        fa = fd
      d = 0.5*(a+b)
      count = count +1
#      print('abs(d-a) = ', abs(d-a))
      
    astar = d
    ier = 0
    return [astar, ier, num]

#Fixed Point Method
def fixedpt(f,x0,tol,Nmax):

    ''' x0 = initial guess''' 
    ''' Nmax = max number of iterations'''
    ''' tol = stopping tolerance'''

    count = 0
    while (count <Nmax):
       count = count +1
       x1 = f(x0)
       if (abs(x1-x0) <tol):
          xstar = x1
          ier = 0
          return [xstar,ier]
       x0 = x1

    xstar = x1
    ier = 1
    return [xstar, ier]

#Secant Method
def secant_method(f,x0,x1,tol,nmax,verb=False):
    #secant (quasi-newton) method to find root of f starting with guesses x0 and x1

    #Initialize iterates and iterate list
    xnm=x0; xn=x1;
    rn=np.array([x1]);
    # function evaluations
    fn=f(xn); fnm=f(xnm);
    msec = (fn-fnm)/(xn-xnm);
    nfun=2; #evaluation counter nfun
    dtol=1e-10; #tolerance for derivative (being near 0)

    if np.abs(msec)<dtol:
        #If slope of secant is too small, secant will fail. Error message is
        #displayed and code terminates.
        if verb:
            print('\n slope of secant at initial guess is near 0, try different x0,x1 \n');
    else:
        n=0;
        if verb:
            print("\n|--n--|----xn----|---|f(xn)|---|---|msec|---|");

        #Iteration runs until f(xn) is small enough or nmax iterations are computed.

        while n<=nmax:
            if verb:
                print("|--%d--|%1.8f|%1.8f|%1.8f|" %(n,xn,np.abs(fn),np.abs(msec)));

            pn = - fn/msec; #Secant step
            if np.abs(pn)<tol or np.abs(fn)<2e-15:
                break;

            #Update guess adding Newton step, update xn-1
            xnm = xn; #xn-1 is now xn
            xn = xn + pn; #xn is now xn+pn

            # Update info and loop
            n+=1;
            rn=np.append(rn,xn);
            fnm = fn; #Note we can re-use this function evaluation
            fn=f(xn); #So, only one extra evaluation is needed per iteration
            msec = (fn-fnm)/(xn-xnm); # New slope of secant line
            nfun+=1;

        r=xn;

        if n>=nmax:
            print("Secant method failed to converge, niter=%d, nfun=%d, f(r)=%1.1e\n'" %(n,nfun,np.abs(fn)));
        else:
            print("Secant method converged succesfully, niter=%d, nfun=%d, f(r)=%1.1e" %(n,nfun,np.abs(fn)));

    return (r,rn,nfun)

#Newton Method
def newton_method(f,df,x0,tol,nmax,verb=False):
    #newton method to find root of f starting at guess x0

    #Initialize iterates and iterate list
    xn=x0;
    rn=np.array([x0]);
    # function evaluations
    fn=f(xn); dfn=df(xn);
    nfun=2; #evaluation counter nfun
    dtol=1e-10; #tolerance for derivative (being near 0)

    if abs(dfn)<dtol:
        #If derivative is too small, Newton will fail. Error message is
        #displayed and code terminates.
        #print('\n derivative at initial guess is near 0, try different x0 \n');
        r = "Broken due to derivative near 0"
        nfun = "None"
    else:
        n=0;
        if verb:
            print("\n|--n--|----xn----|---|f(xn)|---|---|f'(xn)|---|");

        #Iteration runs until f(xn) is small enough or nmax iterations are computed.

        while n<=nmax:
            if verb:
                print("|--%d--|%1.8f|%1.8f|%1.8f|" %(n,xn,np.abs(fn),np.abs(dfn)));

            pn = - fn/dfn; #Newton step
            if np.abs(pn)<tol or np.abs(fn)<2e-15:
                break;

            #Update guess adding Newton step
            xn = xn + pn;

            # Update info and loop
            n+=1;
            rn=np.append(rn,xn);
            dfn=df(xn);
            fn=f(xn);
            nfun+=2;

        r=xn;

        if n>=nmax:
            print("Newton method failed to converge, niter=%d, nfun=%d, f(r)=%1.1e\n'" %(n,nfun,np.abs(fn)));
        else:
            print("Newton method converged succesfully, niter=%d, nfun=%d, f(r)=%1.1e" %(n,nfun,np.abs(fn)));

    return (r,rn,nfun)


#-------Written Problems-------#
def Problem1():
    Ti = 20
    Ts = -15
    alpha = 0.138e-6
    t = 60 * 24 * 3600

    T = lambda x: Ts + (Ti-Ts)*sp.erf(x/(2*np.sqrt(alpha*t)))
    dT = lambda x: (Ti-Ts)*(1/(np.sqrt(np.pi)*np.sqrt(alpha*t)))*np.exp(-((x/(2*np.sqrt(alpha*t)))**2))


    #Part a) Plot, validate xbar
    xbar = 10
    xVec = np.linspace(0,xbar,100)
    yVec = T(xVec)

    plt.figure()
    plt.plot(xVec,yVec)
    plt.grid(True)
    plt.xlabel("X Points [m]")
    plt.ylabel("Temperatures [degrees C]")
    plt.title("1.a - Initial Plot of Temperature Function")
    plt.show()

    #Part b) Bisection Method
    tol = 10**-13
    a = 0
    b = xbar
    [biDepth,_,_] = bisection(T,a,b,tol)
    print("The depth approximated by Bisection was: ", biDepth)

    #Part c) Newton's Method
    nmax = 25
    x0 = 0.01
    [NDepthOG,_,_] = newton_method(T,dT,x0,tol,nmax)
    print("The depth approximated by Newton (x0 = 0.01) was: ", NDepthOG)
    [NDepthBar,_,_] = newton_method(T,dT,xbar,tol,nmax)
    print("The depth approximated by Newton (x0 = x_bar) was: ", NDepthBar)


def Problem5():
    f = lambda x: x**6 - x - 1
    df = lambda x: 6*(x**5)-1
    x0 = 2
    x1 = 1
    tol = 1e-9
    nmax = 20

    [_, rnSecant, _] = secant_method(f,x0,x1,tol,nmax)
    [_, rnNewton, _] = newton_method(f,df,x0,tol,nmax)
    
    bettertol = 1e-16
    [alpha, _, _] = newton_method(f,df,x0,bettertol,nmax)
    
    xVecSecant = []
    yVecSecant = []
    xVecNewton = []
    yVecNewton = []
    for j in range(len(rnSecant)):
        print("Iteration ", j+1, " - Secant Error: ", abs(rnSecant[j]-alpha), " Newton Error: ", abs(rnNewton[j] - alpha))
    
    for k in range(len(rnSecant)-1):
        if k == 1:
            xVecSecant = abs(rnSecant[k] - alpha)
            yVecSecant = abs(rnSecant[k+1] - alpha)
            xVecNewton = abs(rnNewton[k] - alpha)
            yVecNewton = abs(rnNewton[k+1] - alpha)
        else:
            xVecSecant = np.append(xVecSecant,abs(rnSecant[k] - alpha))
            yVecSecant = np.append(yVecSecant,abs(rnSecant[k+1] - alpha))
            xVecNewton = np.append(xVecNewton,abs(rnNewton[k] - alpha))
            yVecNewton = np.append(yVecNewton,abs(rnNewton[k+1] - alpha))

    plt.figure()
    plt.loglog(xVecSecant,yVecSecant)
    plt.loglog(xVecNewton,yVecNewton)
    plt.grid(True)
    plt.xlabel("k+1 Iteration Error")
    plt.ylabel("k Iteration Error")
    plt.legend(["Secant Method", "Newton Method"])
    plt.title("5.b - Convergance Rate Comparison")
    plt.show()

    mSecantAvg = (np.log(yVecSecant[5]) - np.log(yVecSecant[0])) / (np.log(xVecSecant[5]) - np.log(xVecSecant[0]))
    mNewtonAvg = (np.log(yVecNewton[5]) - np.log(yVecNewton[0])) / (np.log(xVecNewton[5]) - np.log(xVecSecant[0]))

    print("The slope of the Secant line is ", mSecantAvg)
    print("The slope of the Newton line is ", mNewtonAvg)

#-------Calling Problem Functions-------#
#Problem 1
Problem1()

#Problem 5
#Problem5()