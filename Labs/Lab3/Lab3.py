import numpy as np
import matplotlib.pyplot as plt

#------------------------------------#
#Given Functions#
#Bisection Method
def bisection(f,a,b,tol):
    
#    Inputs:
#     f,a,b       - function and endpoints of initial interval
#      tol  - bisection stops when interval length < tol

#    Returns:
#      astar - approximation of root
#      ier   - error message
#            - ier = 1 => Failed
#            - ier = 0 == success

#     first verify there is a root we can find in the interval 

    fa = f(a)
    fb = f(b);
    if (fa*fb>0):
       ier = 1
       astar = a
       return [astar, ier]

#   verify end points are not a root 
    if (fa == 0):
      astar = a
      ier =0
      return [astar, ier]

    if (fb ==0):
      astar = b
      ier = 0
      return [astar, ier]

    count = 0
    d = 0.5*(a+b)
    while (abs(d-a)> tol):
      fd = f(d)
      if (fd ==0):
        astar = d
        ier = 0
        return [astar, ier]
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
    return [astar, ier]

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

#------------------------------------#
#Main Code#
#Section 4.1
f = lambda x: (x**2)*(x-1)
tol = 1e-7
a = np.array([0.5,-1,-1])
b = np.array([2,0.5,2])
parts = ["a","b","c"]

#print(a)
#print(b)

for i in range(3):
    [astar, ier] = bisection(f,a[i],b[i],tol)
    print("The root found was: ", astar)
    if (ier == 1):
        print("The root finding method failed for part ", parts[i])
    else:
        print("The root finding method succeded for part ", parts[i])

x_test = np.linspace(-3,3,100)
y_test = f(x_test)

plt.figure()
plt.grid(True)
plt.plot(x_test,y_test)
plt.title("Graphing Plot for part 1")
plt.show()

#Section 4.2
tol = 10^(-5)
#Part a
f_a = lambda x: (x-1)*(x-3)*(x-5)
a_a = 0
b_a = 2.4

[astar,ier] = bisection(f_a,a_a,b_a,tol)
print("The expected zero in this region is 1 and the found zero is ", astar)

#Part b
f_b = lambda x: ((x-1)**2)*(x-3)
a_b = 0
b_b = 2

[astar,ier] = bisection(f_b,a_b,b_b,tol)
print("The expected zero in this region is 1 and the found zero is ", astar)
print("No sign change occured and thus we have an error")

#Part c
f_c = lambda x: np.sin(x)
a_c = 0
b_c = 0.1

[astar,ier] = bisection(f_c,a_c,b_c,tol)
print("The expected zero in this region is 0 and the found zero is ", astar)

a_c = 0.5
b_c = (3*np.pi)/4

[astar,ier] = bisection(f_c,a_c,b_c,tol)
print("There is no expected zero in the region and the found zero is ", astar)

#Section 4.3
g_a = lambda x:  x*((1 + ((7-x**5)/(x**2))) ** 3)
g_b = lambda x: x - ((x**5-7)/(x**2))
g_c = lambda x: x - ((x**5-7)/(5*(x**4)))
g_d = lambda x: x - ((x**5-7)/(12))

#Verification
x_verif = 7**(1/5)
print("Using the given verification x, for each function they return: ",
       [g_a(x_verif),g_b(x_verif),g_c(x_verif),g_d(x_verif)])

#Fixed Point:
#(ran out of time)