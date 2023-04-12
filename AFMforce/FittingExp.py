#!/usr/bin/env python
""" Fitting function to exponential decay
"""

from numpy import abs, sqrt, exp, log, ones, polyfit, polyval, diag
from scipy.optimize import leastsq
from matplotlib import pyplot as pl

def f_exp(x, a, m, y0):
    """ Calculate an exponential decay in the form of:
        a*exp(m*(x)) + y0
        Please consider, that a*exp(m*(x-x0)) + y0 is
        a1 * exp(m*x)+y0, where a1 = a*exp(-m*x0)

        parameters:
        x   an array of x values
        a, m, y0 as in the equation

        return:
        array of y values
        """

    return a*exp(m*x) + y0
#end of f_exp

def J_exp(parms, *args):
    """ Calculate the Jacobian for the leastsq() fitting.

        parms:  x-values
        args:   a, m, y0

        return:
        an array of [dy/da, dy/dm, dy/dy0]
    """
    a, m, F0 = parms
    x = args[0]
    #start calculating:
    dyda = exp(m*x)

    return [ dyda, a*x*dyda, ones(x.shape, dtype='float')]
#end J_exp

def err_exp(parms, *args):
    """ Calculate the chi error for the estimated parameters
        parms: a, m, x0, F0
        args: x,y values, weight values

        return:
            (fitted - y)*weight
    """
    NA = len(args)
    N = len(parms)

    if NA < 2 or NA > 3:
        raise ValueError("Invalid number of parameters")

    x = args[0]
    y = args[1]

    if NA == 2:
        noweight = True
    elif NA == 3:
        weight = args[2]
        noweight = False
    #end if

    #if the x,y length differs, let it crash
    a, m, y0 = parms
    yy = a*exp(m*x) + y0
    chi = (yy-y) if noweight else (yy-y)*weight

    return chi
#end err_power

def fit_exp(x, y, Fmax = None, withweight= False, verbose= False):
    """ Fit an exponential function to the data.
        Try estimating all parameters automatically.

        x,y:    data set, indentation and force values for example
        Fmax:   up to this force do the fitting
        withweight: use a weighting if set
        verbose:    give graphical feedback

        return:
        a dict structure containing the results
    """

    N = len(x)

    if len(y) != N:
        raise ValueError("Length of x and y do not match")
    #end if mismatch

    #set an index to which data we keep for fitting:
    indx = ones(y.shape, dtype='bool') if Fmax is None else y < Fmax

    if indx.sum() < 3:
        raise ValueError("Too short dataset, N < 3")
    #end if too short N

    xx = x[indx]
    yy = y[indx]
    # print("we had {0}".format(N), end=' ')
    N = len(xx)
    # print("we have new length {0}".format(N))

    F0 = yy.min()
    # we cannot take F0 as is, because then we have some 0 values
    if F0 < 0:
        F0 = 1.1*F0
    else:
        F0 = 0.9*F0

    ly = log((yy-F0))
    lmfit = polyfit(xx,ly,1)
    # print("fitted semi-log data")
    # print(lmfit)

    # pl.clf(); pl.plot(xx,ly,'bo')
    # pl.plot(xx, polyval(lmfit, xx), 'r-')
    # pl.draw()

    # our first estimate:
    m = lmfit[0]
    a = exp(lmfit[1])

    if withweight :
        w = abs(yy)
        fit = leastsq( err_exp, [a, m, F0], (xx,yy, w),\
            Dfun = J_exp,
            ftol= 1E-8,\
            col_deriv= 1, full_output= True)

    else:
        fit = leastsq( err_exp, [a, m, F0], (xx,yy),\
            Dfun = J_exp,
            ftol= 1E-8,\
            col_deriv= 1, full_output= True)

    # print("fit complete")
    # print(fit[0])

    a,m,F0 = fit[0]
    chi2 = (fit[2]['fvec']**2).sum()
    chi2red = chi2/(len(xx) - len(fit[0]))
    msg = fit[3]
    r2 = 1.0 - (chi2 / ((yy - yy.mean())**2).sum())

    if fit[1] is not None:
        cov = fit[1]*chi2red
        da, dm, dF0 = sqrt(diag(cov))
        singular = False
    else:
        da, dm, dF0 = (-1, -1, -1)
        singular = True
        cov = None

    if verbose:
        pl.clf()
        pl.plot(x,y,'bo')
        pl.xlabel('distance, $\mu$m')
        pl.ylabel('force, nN')
        pl.plot(x, f_exp(x, a, m, F0), 'r-')

    return {'x': xx, 'y': yy, 'fitted': f_exp(x, a, m, F0),
            'a': a, 'da': da,
            'm': m, 'dm': dm,
            'F0': F0, 'dF0': dF0,
            'r2': r2, 'chi2': chi2,
            'singular': singular,
            'message': msg, 'cov': cov}
#end fit_power
