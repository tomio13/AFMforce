#!/usr/bin/env python
""" A simple linear fit model for general use. The idea is very basic,
    I am astonished that it is not included in numpy.
    """

from numpy import matrix, asarray, sqrt, arange

from numpy.linalg import pinv
from matplotlib import pyplot as pl

def lm( y, xs, verbose=True):
    """ Linear model fitting (without weights)
        Use this function to fit generalized linear models.

        Parameters:
        y           the vector of desired values, that is the
                    data to pass to
        xs          a list of vectosrrs to fit
                    y = b[0] x[:,0] + b[1] x[:,1] + ...
                    will be fitted.
        return:
        a dict containing:
        'b'         the weight vector
        'res'       residuals
        'r2'        the estimated r^2
        'err'       estimated error of b
        'var'       variance calculated from the inverse matrix
        'Mi'        the inverse matrix

        Note:   there is a 'full fetched' solution as linalg.lstsq,
        which does not provide the inverse matrix or estimation
        of the covariance matrix.

    """
    Np = len(xs)
    N = len(y)

    #generate the matrix:
    M = matrix(xs).T
    #the generalized problem solves the (A.T A)b = A.T y equation
    #this A.T A step makes the problem matrix much easier:
    Mi = pinv(M.T*M)
    #for weights it would be: Mi = linalg.pinv(M.T*W*M)

    #make y usable for matrix op:
    y.shape = [y.shape[0],1]

    #with weights it would be: beta = Mi*M.T*W*y
    #the solution:
    beta = (Mi*M.T*y)

    #since the original equation: M beta = y, our fitted y:
    yfitted = asarray( M*beta)
    res = (yfitted - y)**2
    chi2 = res.sum()

    r2 = 1.0 - res.sum()/((y-y.mean())**2).sum()
    #the covariance:
    cov = chi2/(N- Np) *Mi
    varb = cov.diagonal()
    varb.shape = [varb.size,]
    errb = sqrt(varb)

    if verbose:
        print("Fitted number of variables: %d" %len(xs))
        print("Residual error: %.4f" %res.sum())
        print("Regression coefficient R^2: %.4f" %r2)
        print("Results:", beta)
        print("Errors:", errb)
        print("variances:", varb)

        x = arange(len(y))
        pl.clf()
        pl.plot(x,y,'go',alpha=0.3)
        pl.plot(x, yfitted, 'r-')
        pl.plot(x, y-yfitted,'bo',alpha=0.3)
    #end if

    return {'b': asarray(beta.T)[0], 'fitted':yfitted, 'res': res, \
            'chi2': chi2, 'r2': r2, 'Np': len(xs), \
            'Nres': len(y)-len(xs),\
            'cov': cov, 'var':asarray(varb), 'err': errb}
#end of lm

