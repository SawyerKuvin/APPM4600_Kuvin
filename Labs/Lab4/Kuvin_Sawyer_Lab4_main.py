#The following code is the functions and main script used for APPM4600: Lab 4
#By: Sawyer Kuvin
#Last Edit: 2/4/25 - in lab
#Sections completed: ~2.1, 2.2, 3.1, ~3.2 (~: Part of)

import numpy as np
import matplotlib.pyplot as plt

#----------------------------------------#
#Functions:
#Fixed Point Method (Given)
def fixedpt(f,x0,tol,Nmax):

    ''' x0 = initial guess''' 
    ''' Nmax = max number of iterations'''
    ''' tol = stopping tolerance'''

    count = 0
    x_tests = x0
    while (count <Nmax):
       count = count +1
       x1 = f(x0)
       x_tests = np.append(x_tests, x1) 
       if (abs(x1-x0) <tol):
          xstar = x1
          ier = 0
          return [x_tests,ier]
       x0 = x1

    xstar = x1
    ier = 1
    return [x_tests, ier]

#Aitkens Method Function
def Aitkens (f,x0,tol,Nmax):

    ''' x0 = initial guess''' 
    ''' Nmax = max number of iterations'''
    ''' tol = stopping tolerance'''

    count = 0
    x_tests = x0
    while (count <Nmax):
       count = count +1
       x1 = f(x0)
       x_tests = np.append(x_tests, x1) 
       if (abs(x1-x0) <tol):
          xstar = x1
          ier = 0
       x0 = x1

    xstar = x1
    ier = 1

    x_hats = []
    for i in range(len(x_tests)-3):
        n = i
        nPlus1 = i+1
        nPlus2 = i+2
        
        pn = x_tests[n]
        pnPlus1 = x_tests[nPlus1]
        pnPlus2 = x_tests[nPlus2]
        
        x_hat = pn - (((pnPlus1 - pn)**2)/(pnPlus2-2*pnPlus1+pn))
        x_hats = np.append(x_hats,x_hat)
    
        if (abs(x_hat-x0) <tol):
          return [x_tests, x_hats]
        x0 = x1



#---------------------------------------#
#Main Script
# Part 2.1 and 2.2
f_x = lambda x: (1/2)*x+4
x0 = 1.25
tol = 1e-8
Nmax = 100
[x_vec,_] = fixedpt(f_x,x0,tol,Nmax)

#print(x_vec)

g_x = lambda x: (10/(x+4))**(1/2)
p0 = 1.5
p = 1.3652300134140976
tol = 1e-10
Nmax = 100
[p_vec,_] = fixedpt(g_x,p0,tol,Nmax)
print(p_vec)

#finding alpha
alphas = []
for i in range(len(p_vec)-1):
    if i >= 2: 
        n = i
        nPlus1 = i+1
        nMinus1 = i-1
        nMinus2 = i-2
        alphaNum = np.log10(abs(p_vec[nPlus1]-p_vec[n])/abs(p_vec[n]-p_vec[nMinus1]))
        alphaDem = np.log10(abs(p_vec[n]-p_vec[nMinus1])/abs(p_vec[nMinus1]-p_vec[nMinus2]))
        alpha = alphaNum/alphaDem
        
        alphas = np.append(alphas,alpha)

print(alphas)

#finding constant - failed
'''lambs = []
for i in range(len(p_vec)-1):
    if i >= 11: 
        n = i
        nPlus1 = i+1
        lambNum = abs(p_vec[nPlus1]-p)
        lambDem = abs(p[n]-p)
        lamb = lambNum/lambDem
        
        lambs = np.append(lambs,lamb)

print(lambs)'''
        


# Part 3.1
#Testing function:

#Aitken's Method
g_x = lambda x: (10/(x+4))**(1/2)
p0 = 1.5
tol = 1e-10
Nmax = 100
[x_vec,_] = fixedpt(g_x,p0,tol,Nmax)
print(x_vec)
[x_vec,x_hat] = Aitkens(g_x,p0,tol,Nmax)
print(x_hat)


