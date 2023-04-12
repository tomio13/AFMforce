#!/usr/bin/env python
""" Fitting function to generalized power law indentation
"""

from numpy import abs, sqrt, exp, log, ones, polyfit, polyval, diag
from scipy.optimize import leastsq
from matplotlib import pyplot as pl

def f_power(x, a, m, x0, F0):
    """ Calculate a general power law in the form of:
        a*(x-x0)**m + F0

        parameters:
        x   an array of x values
        x0, a, m, F0 as in the equation

        return:
        array of y values
        """

    x = x-x0

    return( a*(x*(x>0))**m + F0)
#end of f_power

def J_power(parms, *args):
    """ Calculate the Jacobian for the leastsq() fitting.

        parms:  x-values
        args:   a, m, x0, F0

        return:
        an array of [dy/da, dy/dm, dy/dx0, dy/dF0]
    """
    a, m, x0, F0 = parms
    x = args[0]
    #start calculating:
    x = x - x0
    dx = (x > 0) if not (x < 0).all() and m >= 1 else 1000*ones(x.shape, dtype='float')

    #dy/da = (x-x0)**m
    dyda = (x*dx)**m
    #dy/dm = a*log(a)*(x-x0)**m
    dydx = -a*m*(x*dx)**(m-1)

    return [ dyda, a*log(dx*x)*dyda, dydx, ones(x.shape, dtype='float')]
#end J_power

def err_power(parms, *args):
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
    a, m, x0, F0 = parms
    x = x - x0
    #only positive parts have a meaning and being evaluated:
    x = x*(x>0)
    if x.sum() == 0:
        return 1000*ones(x.shape, dtype='float')

    yy = a*x**m + F0
    chi = (yy-y) if noweight else (yy-y)*weight

    return chi
#end err_power

def fit_power(x, y, Fmax = None, withweight= False, verbose= True):
    """ Fit a general power law to the data where x > x0.
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
    print("we had {0}".format(N), end=' ')
    N = len(xx)
    print("we have new length {0}".format(N))

    F0 = yy.min()
    #if we have some x0 in the x-axis, we have
    #positive and negative values.
    #if we have only negative, that is a problem,
    #then try correcting it:
    i0 = ( yy == yy.min()).nonzero()[0]
    if len(i0) < 1:
        i0 = 0
    else:
        i0 = i0[0]

    x0 = xx[i0] if xx.all() < 0 else 0
    #we start a log-log fit on the upper slope:
    #indxfit = (xx > x0) & (yy > 0.5*yy.max())
    print("x0:", x0)

    fitindx = (xx >0) & (yy > F0)
    #ly = log((yy-F0)[indxfit])
    #lx = log(xx[indxfit] )
    ly = log( (yy-F0)[fitindx])
    lx = log( xx[fitindx] )
    print("We fit log-log to {0} data".format(len(ly)))

    lmfit = polyfit(lx,ly,1)
    print("fitted log-log data")
    print(lmfit)

    pl.clf(); pl.plot(lx,ly,'bo')
    pl.plot(lx, polyval(lmfit, lx), 'r-')
    pl.draw()

    #our first estimate:
    m = lmfit[0]
    a = exp(lmfit[1])

    if withweight :
        w = abs(yy)* (xx > 0)
        fit = leastsq( err_power, [a, m, x0, F0], (xx,yy, w),\
            Dfun = J_power, ftol= 1E-8,\
            col_deriv= 1, full_output= True)

    else:
        fit = leastsq( err_power, [a, m, x0, F0], (xx,yy),\
            Dfun = J_power, ftol= 1E-8,\
            col_deriv= 1, full_output= True)

    print("fit complete")
    print(fit[0])

    a,m,x0,F0 = fit[0]
    chi2 = (fit[2]['fvec']**2).sum()
    chi2red = chi2/(len(xx) - len(fit[0]))
    msg = fit[3]
    r2 = 1.0 - (chi2 / ((yy - yy.mean())**2).sum())

    if fit[1] is not None:
        cov = fit[1]*chi2red
        da, dm, dx0, dF0 = sqrt(diag(cov))
        singular = False
    else:
        da, dm, dx0, dF0 = (-1, -1, -1, -1)
        singular = True
        cov = None

    pl.figure(2)
    pl.clf()
    pl.plot(x,y,'bo')
    pl.plot(x, f_power(x, a, m, x0, F0), 'r-')

    return {'x': xx, 'y': yy, 'fitted': f_power(x, a, m, x0, F0),
            'a': a, 'da': da,
            'm': m, 'dm': dm,
            'x0': x0, 'dx0': dx0,
            'F0': F0, 'dF0': dF0,
            'r2': r2, 'chi2': chi2,
            'singular': singular,
            'message': msg, 'cov': cov}
#end fit_power
