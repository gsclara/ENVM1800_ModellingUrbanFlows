from pprint import pprint
import numpy as np

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

