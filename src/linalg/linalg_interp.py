import numpy as np

def gauss_iter_solve(A,b,x0=None,tol= 1e-8,alg= 'seidel'):
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

    A = np.array(A,dtype=float)
    b = np.array(b,dtype=float)

    if not (n := A.shape[0]) == (m := A.shape[1]):
        raise ValueError("A's rows {n} and columns {m} are not equal in size")
    if not (dimA := len(A.shape)) == 2:
        raise ValueError("A's is {dimA} dimension but needs to be 2D")
    if not (dimb := len(b.shape)) == (2 or 1):
        raise ValueError("b's is {dimb} dimension but needs to be 1D or 2D")

    n_b = b.shape[0]
    if rs_1d := (dimb == 1):
        b = np.reshape(b, (nb,mb :=1))
    else: mb=b.shape[1]

    if(mb := b.shape[0]) != n:
        raise ValueError("b has {mb} # of rows which should be the same # of rows from A")

    if x0 is None:
        x0 = np.zeros(n).reshape(-1,1)
    if len(x0.shape)==1:
        x0 = np.reshape(x0,(n,1))
    if not rs_1d and x0.shape[1]==1:
        x0 = x0 @ np.ones((1,mb))

    alg= alg.strip().lower()

    eps_a = 1
    eps_s = tol
    i = 1
    max_i = 100

    Adiag = np.diag(np.diag(A))
    Adiaginv = np.linalg.inv(Adiag)
    A0 = A - Adiag
    A0st = Adiaginv @ A0
    bst = Adiaginv @ b

    if alg == 'jacobi':
        while np.max(eps_a) > eps_s and i < max_i:
            x_last = np.array(x0)
            x0 = bst - (A0st @ x0)
            i += 1
            delta_x = x0 - x_last
            eps_a = np.linalg.norm(delta_x, axis =0) / np.linalg.norm(x0, axis = 0)
            if i > max_i:
                raise RuntimeWarning('system not converging')
    elif alg == 'seidel':
        while np.max(eps_a) > eps_s and i < max_i:
           x_last = np.array(x0)
           for i, _ in enumerate(A):
               x0[i,:] = bst[i:i+1,:] - A0st[i:i+1,:] @ x0[:,:]
           i += 1
           delta_x = x0 - x_last
           eps_a = np.linalg.norm(delta_x, axis =0) / np.linalg.norm(x0, axis = 0)
           if i > max_i:
               raise RuntimeWarning('system not converging')

    else:
        raise ValueError('alg must be either jacobi or seidel')
    return x0

def spline_function(xd,yd,order = 3):
    """
    Generates a spline function given two vectors x and y

    Parameters
    -------------------------------------------------------
    xd : array like float [data with increasing value]

    yd : array_like float []

    order : int [default 3 can be 1,2, or 3]

    Returns
    -------------------------------------------------------
    function that take one parameter which is a float or array like float and returns interpolated y value(s)

    Raises
    -------------------------------------------------------
    Value Errors:
        # of independent variable and independ variable are not the same
        xd array contains repeated values
        xd values are not in increasing order

    """
    xd = np.array(xd,dtype=float)
    yd = np.array(yd,dtype=float)

    if len(xd) != len(yd):
        raise ValueError('the # of independent variable and independ variable are not the same')
    xd_unique = np.unique(xd)
    if len(xd_unique) != len(xd):
        raise ValueError('xd array contains repeated values')
    if (all (xd[i] <= xd[i-1] for i in range (len(xd)-1))):
        raise ValueError('xd values are not in increasing order')

    dx = np.diff(xd)
    dy = np.diff(yd)
    dy_dx = dy/dx

    if order == 1:
        def s1(x):
            a = yd[:-1]
            b = dy_dx[:-1]
            for n in x:
                k = np.array([np.nonzero(xd<=n)[0][0]-1 for n in x])
                k = np.where(k<0,0,k)
                y = a[k-1] + b[k-1]*(x-xd[k-1])
            return y
        return s1
    elif order == 2:
        def s2(x):
            num = xd.shape[0]
            right = np.zeros(num-1)
            right[1:] = np.diff(dy_dx, axis=0)
            num = len(xd)
            A = np.zeros((num-1,num-1))
            A[0,0:2] = [1, -1]
            A[1:,:-1] += np.diag(dx[:-1])
            A[1:,1:] += np.diag(dx[1:])
            c = np.linalg.solve(A,right)
            b = dy - (c * dx)
            a = yd[:-1]
            for m in x:
                k = np.array([np.nonzero(xd>=m)[0][0]-1 for m in x])
                k = np.where(k<0,0,k)
                y = a[k] + b[k] * (x-xd[k]) + c[k] * (x - xd[k])**2
                return y
        return s2
    elif order ==3:
        def s3(x):
            num = xd.shape[0]
            diff_matrix = np.diff(dy_dx)
            right = np.zeros(num)
            right[1:-1] = 3 * diff_matrix
            num = len(xd)
            A = np.zeros((num,num))
            A[1,0] = dx[0]
            A[-2,-1] = dx[-1]
            A[0,:3] = [-dx[1], (dx[0]+dx[1]), -dx[0]]
            A[-1,-3:] = [-dx[-1], (dx[-1]+dx[-2]), -dx[-2]]
            A[1:-1,:-2] += np.diag(dx[:-1])
            A[1:-1,1:-1] +=np.diag(2*(dx[:-1]+dx[1:]))
            A[1:-1,2:] += np.diag(dx[1:])

            c = np.linalg.solve(A,right)
            d = np.diff(c)/(dx *3)
            b = dy_dx - dx * ( c[:-1] + c[1:] * 2)/3
            k = np.array([np.nonzero(xd>=n)[0][0]-1 for n in x])
            k = np.where(k<0,0,k)
            y = np.array([(yd[k]+ b[k] * (n -xd[k])+ c[k] * (n - xd[k])**2+ d[k] * (n - xd[k])**3) for k, n in zip(k, x)])
            return y
        return s3
    else:
        raise ValueError('order is not 1, 2, or 3')
