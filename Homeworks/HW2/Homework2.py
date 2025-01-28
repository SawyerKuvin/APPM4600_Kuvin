import random
import numpy as np
import matplotlib.pyplot as plt

#----------------------------------------#
#Problem Functions

#Problem 1

#Problem 2

#Problem 3

#Problem 4
def Problem4a():
    t = np.arange(0,(31*np.pi/30),(np.pi/30))
    #print(t)
    y = np.cos(t)
    s = t * y
    S = sum(s)
    print("the sum is: ", S)

def Problem4b():
    #Part 1
    theta = np.linspace(0,2*np.pi, 100)
    R = 1.2
    dr = 0.1
    f = 15
    p = 0
    x = R * (1 + dr * np.sin(f*theta + p)) * np.cos(theta)
    y = R * (1 + dr * np.sin(f*theta + p)) * np.sin(theta)

    #Plotting Initial Figure
    plt.figure()
    plt.plot(x,y)
    plt.title("4b - Initial Parametric Curve")
    plt.xlabel("X points")
    plt.ylabel("Y points")
    plt.grid(True)
    plt.show()

    #Part 2
    plt.figure()
    for i in range(10):
        R = i
        dr = 0.05
        f = 2+i
        p = random.uniform(0,2)
        x = R * (1 + dr * np.sin(f*theta + p)) * np.cos(theta)
        y = R * (1 + dr * np.sin(f*theta + p)) * np.sin(theta)
        plt.plot(x,y) #Calling plot within the loop
    plt.title("4b - Looped Parametric Curves")
    plt.xlabel("X points")
    plt.ylabel("Y points")
    plt.xlim([-10,10])
    plt.ylim([-10,10])
    plt.grid(True)
    plt.show()




#----------------------------------------#
# Main Script

#Calling Problem 1

#Calling Problem 2

#Calling Problem 3

#Calling Problem 4
Problem4a()
Problem4b()

