##APPM 4600 - Lab 5 (In Lab Code Attempt)
#Author: Sawyer Kuvin
#Date: 2/18/25
#Questions attempted/completed: 2.1, 2.2, ~3.1
# (~ partially completed/attempted) 

import numpy as np
from numpy import linalg as lg
import matplotlib.pyplot as plt
import scipy.special as sp

disp = False

#-------Calculation functions-------#
def Forward_Difference(f,s,h):
    fprime = (f(s+h)-f(s))/h
    return(fprcime)

def Centered_Difference(f,s,h):
    fprime = (f(s+h)-f(s-h))/(2*h)
    return(fprime)


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


    



#-------Running Lab Questions-------#
BeforeLab()