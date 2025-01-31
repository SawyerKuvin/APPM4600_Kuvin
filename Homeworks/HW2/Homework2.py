import random
import numpy as np
import matplotlib.pyplot as plt

#----------------------------------------#
#Problem Functions

#Problem 1

#Problem 2
def Problem2():
    #Part b
    A = (1.0/2.0) * np.array([[1.0,1.0],
                              [1+10**(-10),1-10**(-10)]])
    b = np.array([[1.0],
                  [1.0]])
    x = np.array([[1.0],
                  [1.0]])
    A_inv = np.array([[1-10**(10),10**10],
                      [1+10**(10),-10**(10)]])

    _,svds,_ = np.linalg.svd(A,full_matrices=False)
    #print(svds)

    sigma_max = max(svds)
    sigma_min = min(svds)
    Condition_Number = sigma_max/sigma_min

    print("The Condition Number for A is ", Condition_Number)

    #Part c
    delta_b = np.array([[10^(-5)],[3*10^(-5)]])
    #delta_b = np.array([[10^(-5)],[10^(-5)]])
    #delta_b = np.array([[10^(-2)],[3*10^(-2)]])
    b_pert = b+delta_b
    x_pert = A_inv@b_pert

    err_abs = np.linalg.norm(x-x_pert)
    err_rel = np.linalg.norm(x-x_pert)/np.linalg.norm(x)

    err_rel_in = np.linalg.norm(delta_b)/np.linalg.norm(b)

    print("The absolute error resultant from the perturbation ", err_abs)
    print("The relative error resultant from the perturbation ", err_rel)
    print("The relative error in the inputs is ", err_rel_in)
    print("The Condition number times thr relative error in the inputs is: ", err_rel_in * Condition_Number)

    '''kappa_denom = (np.sqrt(b_pert[0]**2+b_pert[1]**2)/np.sqrt(b[0]**2+b[1]**2))
    kappa_num = ((np.sqrt(err_abs[0]**2+err_abs[1]**2)/np.sqrt(x[0]**2+x[1]**2)))
    kappa = kappa_num/kappa_denom'''

#Problem 3
def Problem3():
    #Justifying part b)
    x_test = np.linspace(-1*10**(-14),1*10**(-14),100)
    y_test = np.e**x_test
    #print((y_test - 1))

    #part c)
    x = 9.999999995000000*(10**(-10))
    y = np.e**x
    print(y-1)

    #part e)
    y_taylor_1 = x + (x**2)/(1*2) #+ (x**3)/(1*2*3) + (x**4)/(1*2*3*4) + (x**5)/(1*2*3*4*5)

    y_taylor = y_taylor_1
    print(y_taylor)




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

#Calling Problem 2
Problem2()

#Calling Problem 3
Problem3()

#Calling Problem 4
Problem4a()
Problem4b()

