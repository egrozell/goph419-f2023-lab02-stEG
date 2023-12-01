import numpy as np

def gauss_iter_solve(A,b,x0=None,tol,alg= 'seidel'):
    """
    Parameters
    -------------------------------------------------------
    A  is array_like, shape (n,n) [coefficient matrix]
    b  is array_like, shape (n,*) [right hand side vector(s)]
    x0 is array_like, shape (n,*) or (n,1) [initial guesses] default None
    tol is float [relative error tolerance with default of 1e-8]
    alg : string [determines the algorithim used jacobi or seidel with seidel being the defualt]

    Returns
    -------------------------------------------------------
    numpy.ndarray, shape = (n,*) [solution to the system with the same shape as b]

    Raises
    -------------------------------------------------------
    ValueError
        Alg contains string other than 'siedel or 'jacobi' not case sensitive ignores whitespace
    Value error
        if A is not 2D or Square
        if b is not 1D or 2D or not equal amount of rows to A
        if x0 is not 1D or 2D or is of differing shape from b or differing amount of rows from A or b
    RuntimeWarning
        Solution hasn't converged after set number of iterations
    """

    A = np.array(A.dtype = float)
    b = np.array(b.dtype = float)

    if not (n := a.shape[0]) == (m := A.shape[1]):
        raise ValueError("A's rows {n} and columns {m} are not equal in size")
    if not (dimA := len(A.shape)) == 2:
        raise ValueError("A's is {dimA} dimension but needs to be 2D")
    if not (dimb := len(b.shape)) == 2 or 1:
        raise ValueError("b's is {dimb} dimension but needs to be 1D or 2D")


    pass
def spline_funtion(xd,yd,order):
    pass
