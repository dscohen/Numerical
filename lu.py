from pylab import *
from myarray import *

'''
This code was adapted from the Matlab given in chapter 5 of
"A First Course on Numerical Methods," by Uri Ascher and Chen Greif.
'''

def ainvb(A,b, bandwidth):
    '''
    Solves the system of equations Ax = b,
    returning 'x'.
    '''
    # compute A = P^T LU
    p,LU = plu(A, bandwidth)
    # solve   
    y = forsub(LU,b,p)
    x = backsub(LU,y)
    return x

#MCC3165489

def plu(A):
    '''
    Perform LU decomposition with scaled partial pivoting.
    Upon returnt the coefficients of L and U replace those
    of the input n-by-n nonsingular matrix A. The row interchanges
    performed are recording in the vector p.
    '''

    LU = A.copy()
    n = shape(LU)[0]
    
    # find scales, initialize permutation vector p
    s = abs(LU).max(1)
    p = arange(n)
    
    # LU decomposition with partial pivoting
    for k in range(n-1):
        # find row index of relative max in column k
        find = 2*bandwidth+1
        q = argmax( abs(LU[k:k+find,k]) / s[k:k+find] )
        q = q + k
        
        # interchange rows k and q and record in p
        t = LU[k,:].copy()
        ts = s[k]
        tp = p[k]
        LU[k,:] = LU[q,:]
        s[k] = s[q]
        p[k] = p[q]
        LU[q,:] = t
        s[q] = ts
        p[q] = tp

        # compute corresponding column of L
        LU[k+1:k+find,k] = LU[k+1:k+find,k] / LU[k,k]

        # update submatrix by outer product
        LU[k+1:k+find,k+1:k+find] = LU[k+1:k+find,k+1:k+find] - outer( LU[k+1:k+find,k], LU[k,k+1:k+find] )

    return (p,LU)

def forsub(A,b,p):
    '''
    Given a unit lower triangular, nonsingular n-by-n matrix A,
    an n-vector b,
    and a permutation p,
    return vector y which solves Ay = Pb.
    '''
    
    n = shape(A)[0]
    
    # permute b according to p
    y = b.copy()
    perm_b = y[p]
    
    # forward substitution
    for k in range(n):
        #limit range to only bandwidth
        y[k] = perm_b[k] - dot( A[k,:k], y[:k] )
    return y

def backsub(A,b):
    '''
    Given an upper triangular, nonsingular n-by-n matrix A
    and an n-vector b,
    return vector x which solves Ax = b
    '''
    
    n = shape(A)[0]
    x = b.copy()
    for k in range(n)[::-1]:
        x[k] = ( b[k] - dot( A[k,k+1:], x[k+1:] ) ) / A[k,k]
    return x

def getL(LU):
    n = shape(LU)[0]
    L = copy(LU)
    for i in range(n):
        L[i,i] = 1
        L[i,i+1:] = 0
    return L

def getU(LU):
    n = shape(LU)[0]
    U = copy(LU)
    for i in range(n):
        U[i,:i] = 0
    return U

def getP(p):
    P = empty( (len(p),len(p)) )
    I = identity( len(p) )
    for i, val in enumerate( p ):
        P[i,:] = I[val,:]
    return P

def test1():
    print '==== test1() ===='
    
    n = 5
    
    import numpy.linalg
    det = 0
    # make a random matrix
    while det < 1e-5:
        diagonals = [rand(n-2), rand(n-1), rand(n), ran(n-1), rand(n-2)]
        A = sp.diags(diagonals, [-2, -1, 0, 1, 2]).todense()
        det = linalg.det(A)
    
    Ainv = zeros((n,n))
    I = identity(n)
    for i in range(n):
        Ainv[:,i] = ainvb(A,I[i,:])
    
    print 'A:'
    print A.round(5)
    print 'Ainv:'
    print Ainv.round(5)
    print 'A * Ainv:' ## Should be the identity matrix
    AAinv = dot( A, Ainv )
    print AAinv.round( 5 )
    
    ## Should be close to zero
    print 'I - A*Ainv:', ( ( I - AAinv ) ** 2 ).sum().sum()


def test2():
    print '==== test2() ===='
    
    n = 5
    
    import numpy.linalg
    det = 0
    while det < 1e-5:
        A = rand(n,n)
        det = linalg.det(A)
    
    ## NOTE: plu() creates LU in-place in its input parameter.
    ##       We still want to have A, so we copy it when passing it.
    p,LU = plu(A.copy())
    print 'p:', p
    print 'LU:'
    print LU.round(5)
    print 'L:'
    print getL( LU ).round(5)
    print 'U:'
    print getU( LU ).round(5)
    print 'A:'
    print A.round(5)
    print 'P.T * L * U:'
    PtLU = dot( dot( getP(p).T, getL(LU) ), getU(LU) )
    print PtLU.round(5)
    print '| A - ( P.T * L * U ) |^2:', ((PtLU - A) ** 2).sum().sum()**2


def main():
    test1()
    test2()

# if the module is executed rather than imported
# run tests
if __name__ == '__main__': main()
