import numpy as np
import matplotlib.pyplot as plt

#python3 Kuvin_Sawyer_Lab1.py | tee Kuvin_Sawyer_Lab1.txt

print("---Section 3.1---")
x = np.linspace(1,20,num=20)
y = np.arange(10,30,2)
print("x equals = ", x)
print("y equals = ", y)
print("Both arrays have 13 elements")

print("Section 3.2 and 3.3")
x3 = x[0:3]
print("The first three entries of x are ", x3)

print("Section 3.4 and 3.5")
w = 10*(-np.linspace(1,10,10))
x2 = np.linspace(1,10,len(w))
plt.semilogy(x2,w)
plt.xlabel("Integer Values")
plt.ylabel("W entries")
