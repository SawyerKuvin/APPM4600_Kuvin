##APPM 4600 - Lab 5 (In Lab Code Attempt)
#Author: Sawyer Kuvin
#Date: 2/11/25
#Questions attempted/completed: 3.2,~3.6a,~3.6b,~3.6c
# (~ partially completed/attempted) 


import numpy as np
import matplotlib.pyplot as plt

verb = 0;

def BiNewton(f,df,a,b,tolBi,tolNew,Nmax):
    fa = f(a)
    fb = f(b)
    count = 0
    if (fa*fb>0):
       ier = 1
       astar = a
       return [astar, ier,count]

#   verify end points are not a root 
    if (fa == 0):
      astar = a
      ier =0
      return [astar, ier,count]

    if (fb ==0):
      astar = b
      ier = 0
      return [astar, ier,count]

    d = 0.5*(a+b)
    rn=np.array([d]);
    while (abs(d-a)> tolBi):
      fd = f(d)
      if (fd ==0):
        astar = d
        ier = 0
        return [astar, rn,count]
      if (fa*fd<0):
         b = d
      else: 
        a = d
        fa = fd
      d = 0.5*(a+b)
      rn = np.append(rn,d)
      count = count +1
#      print('abs(d-a) = ', abs(d-a))
#newton method to find root of f starting at guess x0
#Initialize iterates and iterate list

##Beginning Newton Method
    xn=d;
    r=d;
# function evaluations
    fn=f(xn); dfn=df(xn);
    nfun=2; #evaluation counter nfun
    dtol=1e-10; #tolerance for derivative (being near 0)
    if abs(dfn)<dtol:
    #If derivative is too small, Newton will fail. Error message is
    #displayed and code terminates.
        print("WE'RE IN NEWTON BABY")
        if verb:
            print('\n derivative at initial guess is near 0, try different x0 \n');

        else:
            n=0;
            if verb:
                print("\n|--n--|----xn----|---|f(xn)|---|---|f'(xn)|---|");
    #Iteration runs until f(xn) is small enough or nmax iterations are computed.
        while n<=Nmax:
            #if verb:
                #print("|--%d--|%1.8f|%1.8f|%1.8f|" %(n,xn,np.abs(fn),np.abs(dfn)));
            
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
        if n>=Nmax:
            print("Newton method failed to converge, niter=%d, nfun=%d, f(r)=%1.1e\n'" %(n,nfun,np.abs(fn)));
        else:
            print("Newton method converged succesfully, niter=%d, nfun=%d, f(r)=%1.1e" %(n,nfun,np.abs(fn)));
            return [r,rn,nfun+count]

#-----------Main Script-------------#

#3.6 Part c)
f = lambda x: np.exp(x**2+7*x-30) - 1
df = lambda x: (2*x+7) * np.exp(x**2+7*x-30)
Nmax = 50
tolBi = 1
tolNew = 1e-10
a = 2
b = 4.5
[r,rn,numb] = BiNewton(f,df,a,b,tolBi,tolNew,Nmax)
print(r,numb)