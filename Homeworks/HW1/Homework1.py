import numpy as np
import matplotlib.pyplot as plt
#The following code is done written within this file and then
#run in the terminal


# PROBLEM 1

x = np.arange(1.920, 2.080,0.001)
p1 = x**9 - 18*(x**8) + 144*(x**7) - 672*(x**6) + 2016*(x**5) - 4032*(x**4) + 5376*(x**3) - 4608*(x**2) + 2304*x - 512
plt.plot(x,p1)
plt.xlabel("X Points Given")
plt.ylabel("P Points")
plt.title("p(x) Plotted in Coefficient Form")
plt.show()

p2 = (x-2)**9
plt.plot(x,p1)
plt.xlabel("X Points Given")
plt.ylabel("P Points")
plt.title("p(x) Plotted in Original Binomial Form")
plt.show()

# PROBLEM 2
x = np.linspace(1.99999999,2.00000001,1000000)
y = np.linspace(0.5,2.5,1000000)
diff = np.abs(x-y)
z1_bad = np.sin(x) - np.sin(y)
z1_good = 2 * np.cos((x+y)/2) * np.sin((x-y)/2)
plt.plot(diff,z1_bad)
plt.plot(diff,z1_good)
plt.show()