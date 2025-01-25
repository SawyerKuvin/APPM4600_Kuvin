import numpy as np
import matplotlib.pyplot as plt
#The following code is done written within this file and then
#run in the terminal


# PROBLEM 1

x = np.arange(1.920, 2.080,0.001)
p1 = x**9 - 18*(x**8) + 144*(x**7) - 672*(x**6) + 2016*(x**5) - 4032*(x**4) + 5376*(x**3) - 4608*(x**2) + 2304*x - 512
plt.figure(1)
plt.plot(x,p1)
plt.xlabel("X Points Given")
plt.ylabel("P Points")
plt.title("p(x) Plotted in Coefficient Form")
plt.show()

p2 = (x-2)**9
plt.figure(2)
plt.plot(x,p2)
plt.xlabel("X Points Given")
plt.ylabel("P Points")
plt.title("p(x) Plotted in Original Binomial Form")
plt.show()

# PROBLEM 2
#Failed attempt to manually determine the error
'''x = np.linspace(1.99999999,2.00000001,1000000)
y = np.linspace(0.5,2.5,1000000)
diff = np.abs(x-y)
z1_bad = np.sin(x) - np.sin(y)
z1_good = 2 * np.cos((x+y)/2) * np.sin((x-y)/2)
plt.figure(3)
plt.plot(diff,z1_bad)
plt.plot(diff,z1_good)
plt.title("Error Determination Attempt - Failed")
plt.show()'''

# PROBLEM 4
# Using python  to find roots to a precision higher than 3 decimal places
# x^2-56x+1=0
a = 1
b = -56
c = 1

descrim = b**2 - 4*a*c
r1 = (-b + np.sqrt(b**2 - 4*a*c)) / (2 * a)
r2 = (-b - np.sqrt(b**2 - 4*a*c)) / (2 * a)

print("r1 is ", r1)
print("r2 is ", r2)


# PROBLEM 5
# Creating data to be turned into a table for part b) and c)
# where the error in the data will be compred between each set


x1 = (3* np.pi)/4
x2 = 10**7
pow_vec = np.arange(16,-1,-1)
del_vec = (1/10)**(pow_vec)

#print(pow_vec)
#print(del_vec)

# Original formula (a)
fx1_a = np.cos(x1+del_vec)-np.cos(x1)
fx2_a = np.cos(x2+del_vec)-np.cos(x2)

# Trig simplification (b)
fx1_b = -2 * np.sin(x1 + (del_vec / 2)) * np.sin(del_vec / 2)
fx2_b = -2 * np.sin(x2 + (del_vec / 2)) * np.sin(del_vec / 2)


# Plot comparing (a) and (b) for x1
plt.figure(4)
plt.grid(True)
plt.semilogx(del_vec,abs(fx1_a-fx1_b))
plt.xlabel("Delta vector")
plt.ylabel("Difference between models")
plt.title("Difference comparison for Original and Trig Model [3pi/4]")
plt.show()

# Plot comparing (a) and (b) for x2
plt.figure(4)
plt.grid(True)
plt.semilogx(del_vec,abs(fx2_a-fx2_b))
plt.xlabel("Delta vector")
plt.ylabel("Difference between models")
plt.title("Difference comparison for Original and Trig Model [10^7]")
plt.show()



# Iterating most accurate ksi for part (c)
# Test for x1, iterating through fractions of delta to determine a ksi value that
# minimizes error 
fx1_c = np.zeros(17)
del_fracs = np.linspace(0,1,100000)
print(del_fracs)
for d in range(17):
    ksi1 = x1+(del_vec[d]*del_fracs)
    fx1_c_test = -(del_vec[d] * np.sin(x1)) - ((del_vec[d]**2/2) * np.cos(ksi1))
    err = abs(fx1_b[d] - fx1_c_test)
    #print(err)
    min_arg1 = np.argmin(err)
    #print(min_arg1)
    fx1_c[d] = -(del_vec[d] * np.sin(x1)) - ((del_vec[d]**2/2) * np.cos(ksi1[min_arg1]))

print(abs(fx1_c-fx1_b))

diffs1 = abs(fx1_c-fx1_b)
plt.figure(5)
plt.grid(True)
plt.semilogx(del_vec,diffs1)
plt.xlabel("Delta vector")
plt.ylabel("Difference between models")
plt.title("Difference comparison for Trig model and Algorithm [3pi/4]")
plt.show()

#Plot with first poor elements removed
plt.figure(6)
plt.grid(True)
plt.semilogx(del_vec[:-3],diffs1[:-3])
plt.xlabel("Delta vector")
plt.ylabel("Difference between models")
plt.title("Difference comparison for Trig model and Algorithm [3pi/4]")
plt.show()

fx2_c = np.zeros(17)
del_fracs = np.linspace(0,1,100000)
print(del_fracs)
for d in range(17):
    ksi2 = x2+(del_vec[d]*del_fracs)
    fx2_c_test = -(del_vec[d] * np.sin(x2)) - ((del_vec[d]**2/2) * np.cos(ksi2))
    err = abs(fx2_b[d] - fx2_c_test)
    #print(err)
    min_arg2 = np.argmin(err)
    #print(min_arg1)
    fx2_c[d] = -(del_vec[d] * np.sin(x2)) - ((del_vec[d]**2/2) * np.cos(ksi2[min_arg2]))

print(abs(fx2_c-fx2_b))

diffs2 = abs(fx2_c-fx2_b)
plt.figure(7)
plt.grid(True)
plt.semilogx(del_vec,diffs2)
plt.xlabel("Delta vector")
plt.ylabel("Difference between models")
plt.title("Difference comparison for Trig model and Algorithm [10^7]")
plt.show()

#Plot with first poor elements removed
plt.figure(8)
plt.grid(True)
plt.semilogx(del_vec[:-2],diffs2[:-2])
plt.xlabel("Delta vector")
plt.ylabel("Difference between models")
plt.title("Difference comparison for Trig model and Algorithm [10^7]")
plt.show()