import numpy as np
import numpy.linalg as la
import math
def driver():
    x = [[1,0],
        [0,1]]

    y = [[2,1],
        [4,7]]



# evaluate the dot product of y and w
    matmult = matrixMult(x,y)

# Using in house function
    matmult_better = np.matmul(x,y)
# print the output
    print("the dot product is : ", matmult)
    print("the dot product computed by numpy is: ", matmult_better)
    return

## --- Question 3 Answers --- 
# When testing both the dot product and the matrix multiplication functions
# in numpy compared to those that was coded in lab, the numpy functions were
# far more time efficient. This is likely due to the optimized/vectorized
# nature of those functions documentation, compared to the non-ideal code that
# was quickly written below and in testDot.py


def matrixMult(x,y):
# Computes the dot product of the n x 1 vectors x and y
    n1 = np.size(x,0)
    m1 = np.size(x,1)
    n2 = np.size(y,0)
    m2 = np.size(y,1)

    if m1 == n2:
        matmult = np.zeros((n1,m2))
        for i in range(len(x)):
            for j in range(len(y[0])):
                for k in range(len(y)):
                    matmult[i][j] += x[i][k] * y[k][j]        
    else:
        print(":(")

    return matmult
driver()