##APPM 4600 - Homework 10
#Author: Sawyer Kuvin
#Date: 4/4/25

#-------Importing Libraries-------#
import numpy as np
from numpy import linalg as lg
import matplotlib.pyplot as plt
import scipy.special as sp
from scipy.linalg import lu_factor, lu_solve, inv

#-------Problem Functions--------#
def Problem1():
    a = 0
    b = 5
    N = 100
    x_test = np.linspace(a,b,N)
    
    func = lambda x: np.sin(x)
    Mfunc = lambda x: x-(x**3)/6+(x**5)/120
    P3_3func = lambda x: (60*x -7*(x**3))/(60+3*(x**2)) #Note: also the P4_2 Function
    P2_4func = lambda x: (360*x)/(360+60*(x**2)+7*(x**4))

    errM = abs(func(x_test)-Mfunc(x_test))
    errP3_3 = abs(func(x_test)-P3_3func(x_test))
    errP2_4 = abs(func(x_test)-P2_4func(x_test))

    plt.figure()
    plt.plot(x_test,func(x_test),'k-',label="Original Function")
    plt.plot(x_test,Mfunc(x_test),'r--',label='MacClaurin Estimation')
    plt.plot(x_test,P3_3func(x_test),'b--',label='Pade Estimation (parts a and c)')
    plt.plot(x_test,P2_4func(x_test),'g--',label='Pade Estimation (part b)')
    plt.grid(True)
    plt.grid('minor')
    plt.title("Different Estimation methods")
    plt.ylabel("Values of Model and Function")
    plt.xlabel("X points")
    plt.legend()
    plt.show()
    
    
    plt.figure()
    plt.semilogy(x_test,errM,'r--',label='MacClaurin Estimation')
    plt.semilogy(x_test,errP3_3,'b--',label='Pade Estimation (parts a and c)')
    plt.semilogy(x_test,errP2_4,'g--',label='Pade Estimation (part b)')
    plt.grid(True)
    plt.grid('minor')
    plt.title("Error of Different Estimation methods")
    plt.ylabel("Error between Model and Function")
    plt.xlabel("X points")
    plt.legend()
    plt.show()







#-------Calling Problem Functions--------#
Problem1()