from pprint import pprint
import numpy as np
# from numpy import array, zeros, diag, diagflat, dot

def jacobi(A,b,N=25,x=None):
    """Solves the equation Ax=b via the Jacobi iterative method."""
    # Create an initial guess if needed                                                                                                                                                            
    if x is None:
        x = np.zeros(len(A[0]))

    # Create a vector of the diagonal elements of A                                                                                                                                                
    # and subtract them from A                                                                                                                                                                     
    D = np.diag(A)
    R = A - np.diagflat(D)

    # Iterate for N times                                                                                                                                                                          
    for i in range(N):
        x = (b - np.dot(R,x)) / D
    return x

A = np.array([[3.0,-1.0, 0.0, 0.0],[-1.0,2.0,-1.0,0.0],[0.0,-1.0,2.0,-1.0],[0.0,0.0,-1.0,3.0]])
b = np.array([202.5,2.5,2.5,802.5])
guess = np.array([100,100,100,100])

sol = jacobi(A,b,N= 25,x=guess)

print ("A:")
pprint(A)

print ("b:")
pprint(b)

print ("x:")
pprint(sol)

def gauss_seidel(A, b, x0, epsilon, max_iterations):
    n = len(A)
    x = x0.copy()

    #Gauss-Seidal Method [By Bottom Science]

    for i in range(max_iterations):
        x_new = np.zeros(n)
        for j in range(n):
            s1 = np.dot(A[j, :j], x_new[:j])
            s2 = np.dot(A[j, j + 1:], x[j + 1:])
            x_new[j] = (b[j] - s1 - s2) / A[j, j]
        if np.allclose(x, x_new, rtol=epsilon):
            return x_new
        x = x_new
    return x

AA = np.array([[3.0,-1.0, 0.0, 0.0],[-1.0,2.0,-1.0,0.0],[0.0,-1.0,2.0,-1.0],[0.0,0.0,-1.0,3.0]])
bb = np.array([202.5,2.5,2.5,802.5])
x0 = np.array([100,100,100,100])
eps = 1e-5
max_iter = 100

x = gauss_seidel(AA, bb, x0, eps, max_iter)
print(x)