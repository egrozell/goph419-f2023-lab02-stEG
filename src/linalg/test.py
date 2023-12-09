import numpy as np
from linalg import linalg_interp as l

def main():
    parta_data = np.loadtxt('../../data/test.txt')

    x = parta_data[:,0]
    y = parta_data[:,1]

    A = np.array(([[-1, 2,-1, 0, 0, 0, 0, 0, 0, 0, 0],
                   [ 1, 4, 1, 0, 0, 0, 0, 0, 0, 0, 0],
                   [ 0, 1, 4, 1, 0, 0, 0, 0, 0, 0, 0],
                   [ 0, 0, 1, 4, 1, 0, 0, 0, 0, 0, 0],
                   [ 0, 0, 0, 1, 4, 1, 0, 0, 0, 0, 0],
                   [ 0, 0, 0, 0, 1, 4, 1, 0, 0, 0, 0],
                   [ 0, 0, 0, 0, 0, 1, 4, 1, 0, 0, 0],
                   [ 0, 0, 0, 0, 0, 0, 1, 4, 1, 0, 0],
                   [ 0, 0, 0, 0, 0, 0, 0, 1, 4, 1, 0],
                   [ 0, 0, 0, 0, 0, 0, 0, 0, 1, 4, 1],
                   [ 0, 0, 0, 0, 0, 0, 0, 0,-1, 2,-1]]),dtype=float)
    B = np.zeros(11,dtype=float)
    for i in range(1,10):
        a = (y[i-1]-y[i])/(x[i-1]-x[i])
        b = (y[i-2]-y[i-1])/(x[i-2]-x[i-1])
        B[i] = 3*(a-b)
    print(A)
    b = np.resize(B, (11,1))
    print(b)
    xn = l.gauss_iter_solve(A,b,alg='jacobi')
    x1 = np.array(xn[:,1])
    print(x1)
    xb = l.gauss_iter_solve(A,b,alg='seidel')
    x2 = np.array(xb[:,1])
    print(x2)


if __name__ == "__main__":
    main()
