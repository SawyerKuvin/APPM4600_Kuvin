import numpy as np
import numpy.linalg as la
import math
def driver():
    n = 100
    x = np.linspace(0,np.pi,n)
# this is a function handle. You can use it to define
# functions instead of using a subroutine like you
# have to in a true low level language

    # functions listed below have been changed for 4.1 to be orthagonal 
    f = lambda x: 2*np.sin(x) + 1.25
    g = lambda x: 2*np.cos(x)
    y = f(x)
    w = g(x)
# evaluate the dot product of y and w
    dp = dotProduct(y,w,n)
# print the output
    print("the dot product is : ", dp)

#Using in house dot production function 4.3
    dp_better = np.dot(y,w,out=None)
    print("the numpy Dot Product is: ",dp_better)
    return

def dotProduct(x,y,n):
# Computes the dot product of the n x 1 vectors x and y
    dp = 0.
    for j in range(n):
        dp = dp + x[j]*y[j]
    return dp
driver()