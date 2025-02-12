################################################################################
# This python script presents examples regarding the secant method and its
# application to 1D nonlinear root-finding, as presented in class.
# APPM 4650 Fall 2021
################################################################################
# Import libraries
import numpy as np;
import matplotlib.pyplot as plt;

# First, we define a function we will test the secant method with.
# Our test function from previous sections
def fun(x):
    return x + np.cos(x)-3;
def dfun(x):
    return 1 - np.sin(x);

################################################################################
# We now implement the Lazy Newton (Chord) method
def lazynewton_method(f,df,x0,tol,nmax,verb=False):
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
        if verb:
            fprintf('\n derivative at initial guess is near 0, try different x0 \n');
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
            #dfn=df(xn);
            fn=f(xn);
            nfun+=1;

        r=xn;

        if n>=nmax:
            print("Lazy Newton method failed to converge, niter=%d, nfun=%d, f(r)=%1.1e\n'" %(n,nfun,np.abs(fn)));
        else:
            print("Lazy Newton method converged succesfully, niter=%d, nfun=%d, f(r)=%1.1e" %(n,nfun,np.abs(fn)));

    return (r,rn,nfun)
################################################################################
################################################################################
# We now implement the Secant method
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
            fprintf('\n slope of secant at initial guess is near 0, try different x0,x1 \n');
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
################################################################################
################################################################################
# We now implement the Newton method
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
        if verb:
            fprintf('\n derivative at initial guess is near 0, try different x0 \n');
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
################################################################################
# Now, we apply this method to our test function
(rl,rnl,nfunl)=lazynewton_method(fun,dfun,3,1e-14,500,True);
(rNR,rNRn,nfunNR)=newton_method(fun,dfun,3,1e-14,100,True);
(r,rn,nfun)=secant_method(fun,3,4,1e-14,100,True);

print("Secant method log10 errors:")
print(np.log10(np.abs(rn-r+1e-17)));
print("Newton method log10 errors:")
print(np.log10(np.abs(rNRn-r+1e-17)));

# We plot n against log10|f(rn)|
plt.plot(np.arange(0,rnl.size),np.log10(np.abs(rnl-r)),'r-o',label='Lazy Newton');
plt.plot(np.arange(0,rn.size),np.log10(np.abs(rn-r+1e-17)),'k-o',label='Secant');
plt.plot(np.arange(0,rNRn.size),np.log10(np.abs(rNRn-r+1e-17)),'b-o',label='Newton');
plt.xlabel('n'); plt.ylabel('log10|rn-r|');
plt.legend();
plt.suptitle("Comparison of Newton vs Quasi Newton methods");
plt.show();
input();

plt.plot(np.arange(0,rn.size),np.log10(np.abs(rn-r+1e-17)),'k-o',label='Secant');
plt.plot(np.arange(0,rNRn.size),np.log10(np.abs(rNRn-r+1e-17)),'b-o',label='Newton');
plt.xlabel('n'); plt.ylabel('log10|rn-r|');
plt.legend();
plt.suptitle("Comparison of Newton vs Quasi Newton methods");
plt.show();
input();

#plt.plot(np.arange(0,rn.size),np.log10(np.abs(fun(rn))),'r-o');
#plt.xlabel('n'); plt.ylabel('log10|f(rn)|');
#plt.suptitle("Secant method results");
#plt.show();
#input();
################################################################################
# We now test the secant method on the exact same problems Newton struggled.
# A good (but not perfect) rule of thumb is, if Newton struggles, so does secant
# (being a quasi Newton method). Again, it must be paired with safeguards
# to become fully robust.

# Double root example (at x=0)
def fun2(x):
    return np.power(x,2)*(np.cos(x)+2);
def dfun2(x):
    return (2*x)*(np.cos(x)+2) + np.power(x,2)*(-np.sin(x));

(r2,rn2,nfun2)=secant_method(fun2,0.5,1.0,1e-14,100,True);
plt.plot(np.arange(0,rn2.size),np.log10(np.abs(rn2)),'r-o');
plt.xlabel('n'); plt.ylabel('log10|rn|');
plt.suptitle("Secant method results (double root at x=0)");
plt.show();
input();

# Example where Newton has cyclical behavior (cubic)
def fun3(x):
    return np.power(x,3)-2*x+2;
def dfun3(x):
    return 3*np.power(x,2)-2;

(r3,rn3,nfun3)=secant_method(fun3,0.1,1.1,1e-14,100,True);
plt.plot(np.arange(0,rn3.size),np.log10(np.abs(fun3(rn3))),'g-^');
plt.xlabel('n'); plt.ylabel('log10|f(rn)|');
plt.suptitle("Secant method results (Newton was cyclical)");
plt.show();
input();

# Example where Newton took its time to get to basin of quadratic convergence
def fun4(x):
    return x + np.cos(2*x) - 3;
def dfun4(x):
    return 1 - 2*np.sin(2*x);

(r4,rn4,nfun4)=secant_method(fun4,0,0.26,1e-14,200,True);
plt.plot(np.arange(0,rn4.size),np.log10(np.abs(fun4(rn4))),'g-^');
plt.xlabel('n'); plt.ylabel('log10|f(rn)|');
plt.suptitle("Newton method results (slow outside quadratic convergence basin)");
plt.show();
input();

# Example where Newton diverged (derivative has singularity at root)
def fun5(x):
    return np.cbrt(x);
def dfun5(x):
    return (1/3)/np.power(np.cbrt(x),2);

(r5,rn5,nfun5)=secant_method(fun5,0.5,0.25,1e-14,100,True);
plt.plot(np.arange(0,rn5.size),np.log10(np.abs(rn5)),'g-^');
plt.xlabel('n'); plt.ylabel('log10|rn|');
plt.suptitle("Secant method results (Newton diverged from root at x=0)");
plt.show();
input();
