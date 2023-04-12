""" Fit a combined model, using a Hertzian ball model close to the surface and an
    exponential decay (as for simple DLVO) for larger distances
"""

# from AFMforce import *
from AFMforce.AFMforce import Smooth

from numpy import abs, sqrt, pi, tan, ones, zeros, polyfit, polyval, exp, log, diag
from scipy.optimize import leastsq
from matplotlib import pyplot as pl

def f_H_exp(x, A, d0, i0, b, l):
    """ Construct a curve of two functions along x such that from
        0 - i0      y=  A (d0-x)**1.5 and
        i0 - end    y= b exp(-x/l)

        Why to use this simple exponential: if we have a shift
        x0 for the exponential, we get: b exp(-(x-x0)/l) =
        b exp(x0/l) exp(-x/l), thus just b is changed.

        Parameters:
        x       array of x coordinages
                This is the distance from the surface! d0 - x is
                indentation.

        A       the amplitude of Hertz model:
                A = 4/3 sqrt(R) R/(1-\nu**2)
        d0      zero point = apparent contact point for the Hertz model
        i0      a value between 0 and len(x), where the two
                functions swirch role
        b       amplitude of the exponential function
        l       decay length of the exponential function
    """
    N = len(x)
    y = zeros(x.shape)
    if i0 < 0 or i0 >= N:
        return y

    # 1. d0 - x is the indentaiton we need
    # 2. x > d0 is where d0 - x is negative, where y is 0
    i0 = int(i0)
    # x > d0 is the area where indentation is < 0
    # we do not want to reach that point:
    i01 = (x[:i0] > d0).nonzero()[0]
    i01 = i01[0] if i01.sum() > 0 else i0

    if i01 > 0:
        y[:i01] = A * (d0 - x[:i01])**1.5

    y[i0:] = b*exp(-x[i0:]*l)

    return y
# end f_H_exp

def err_H_exp_i0(parms, *args):
    """ calculate the error for the given parameters between the
        provided y array and the f_H_exp() with the parameters.
        In this version, i0 is fixed, not tunable.

        Provide huge errors for cases where:
        l < 0
        b < 0
        A < 0

        params:     A, d0, b, l sent to f_H_exp
        args        i0, x, y, weight

        return:     weight*(f_H_exp(x, A, d0, i0, b, l) - y)

        leastsq will handle the squaring and summing
    """

    NA = len(args)
    Np = len(parms)

    if NA == 3:
        i0 = args[0]
        x = args[1]
        y = args[2]
        noweight = True

    elif NA == 4:
        i0 = args[0]
        x = args[1]
        y = args[2]
        weight = args[3]
        noweight = False
    else:
        raise ValueError("Invalid dataset")


    A, d0, b, l = parms
    N = len(x)

    # penalty for invalid parameter values:
    err = 0
    f = 1E8
    if A < 0:
        err += -f*A
    # this should not happen:
    if i0 < 0:
        err += -f*i0
    elif i0 > N:
        err += f*i0

    if b < 0:
        err += -f*b
    # l is in microns, so in the 10^-3 order
    if l < 0:
        err += -f*1000*l

    if d0 < x[0] or d0 > x[-1]:
        err += abs(d0)*f


    chi = f_H_exp(x, A, d0, i0, b, l) - y

    if noweight:
        return err + chi

    else:
        return err + weight*chi
#end of err_H_exp_i0

def err_H_exp(parms, *args):
    """ calculate the error for the given parameters between the
        provided y array and the f_H_exp() with the parameters.

        Provide huge errors for cases where:
        i0 < 0 or i0 > N-1 (N = len(x))
        l < 0
        b < 0
        A < 0

        params:     A, d0, i0, b, l sent to f_H_exp
        args        x, y, weight

        return:     weight*(f_H_exp(x, A, d0, i0, b, l) - y)

        leastsq will handle the squaring and summing
    """

    NA = len(args)
    Np = len(parms)

    if NA == 2:
        x = args[0]
        y = args[1]
        noweight = True

    elif NA == 3:
        x = args[0]
        y = args[1]
        weight = args[2]
        noweight = False
    else:
        raise ValueError("Invalid dataset")


    A, d0, i0, b, l = parms
    N = len(x)

    # penalty for invalid parameter values:
    err = 0
    f = 1E8
    if A < 0:
        err += -f*A
    if i0 < 0:
        err += -f*i0
    elif i0 > N:
        err += f*i0

    if b < 0:
        err += -f*b
    # l is in microns, so in the 10^-3 order
    if l < 0:
        err += -f*1000*l

    if d0 < x[0] or d0 > x[-1]:
        err += abs(d0)*f


    chi = f_H_exp(x, A, d0, i0, b, l) - y

    if noweight:
        return err + chi

    else:
        return err + weight*chi
# end err_H_exp

def fit_H_exp(x, y,
                R = 1,
                nu= 0.5,
                Flim = 0.5,
                Fnoise = 0.3,
                Rg = 400,
                Wg = 100,
                with_i0 = False,
                verbose= False):
    """ Fit a combined Hertz modela nd exponential function to the data.
        Assume that the data is sorted such the closest points are first.

        The fit assumes that fitting a Hertz model to the highest part of the
        curve (around zero) will have a deviation as the force decreases.
        Find the deviation calculating the absolute error smoothed with a
        Gaussian kernel, and see where this gets higher than the error level
        provided by Fnoise.
        This can be for example 5x standard deviation of the baseline.

        Fit the rest of the curve using an exponential function. While the
        Hertz model tunes an apparent contact point, it is not the same for
        the exponential, because a*exp(-b(x-x0)) = (a*exp(b x0))*exp(-b x).

        parameters
        x, y,:      separation in micrometers and force in nN
        R           radius of the sphere, in micron
        nu          Poisson ratio
        Flim        limiting fraction of the maximal force;
                    the linearized estimate of the Hertz fit is done in the
                    (Flim ... 1)Fmax range
        Fnoise      the noise level of the background in nN
        Rg          window radius for a Gaussian kernel
        Wg          width (sigma) of the Gaussian
        eith_i0     make i0 a fit parameter
        verbose     do some plotting
    """

    if Flim > 1 or Flim <= 0:
        print('Invalid Flim:', Flim, 'resetting to 0.3')
        Flim = 0.3

    # first half is Hertz(?):
    N = len(x)
    Fmax = y.max()
    i1 = (y > Flim*Fmax).nonzero()[0][-1]
    dx = -x[:i1]
    dy = y[:i1]
    indx = dx > 0
    if indx.sum() == 0:
        print('All negative indentation!')
        dx = x[:i1].max() - x[:i1]
    else:
        dx = dx[indx]
        dy = dy[indx]

    dx = dx**1.5
    ft = polyfit(dx, dy, 1)
    A = ft[0]
    d0 = 0
    yfit = zeros(y.shape)
    yfit[ x < 0 ] = polyval(ft, (-x[x < 0])**1.5)

    err = Smooth(abs(y - yfit), 400, 100, 'Gauss')

    indx = (err > Fnoise).nonzero()[0]
    # problems arise if err is too high over all,
    # of err is way too low (Fnoise is too high)
    if len(indx) == 0:
        # all error is below Fnoise, we need a new value
        print('Re-estabilsh i0')
        Fnoise = err[:10].max() if len(err) > 10 else err[0]
        indx = (err > Fnoise).nonzero()[0]
        if len(indx) == 0:
            #now this means the start dominates the error
            # we just pick a point at the end:
            print('Fallback i0')
            i0 = int(0.5*N)
        else:
            i0 = indx[0]
    else:
        i0 = indx[0]

    if verbose:
        pl.clf()
        pl.plot(x, y, '+')
        pl.xlabel('distance, $\mu$m')
        pl.ylabel('force, nN')
        pl.ylim([min(0, y.min()), 1.2*y.max()])
        pl.plot(x, yfit, '-', color= 'orange', linewidth= 2, alpha= 0.5)
        pl.plot(x, err, 'r-', linewidth = 2)
        pl.plot((x.min(),x.max()), (Fnoise, Fnoise), 'r--')

    # validity of the fit when the error goes off the noise
    print('Cutting at:', int(i0), x[i0], y[i0], err[i0])
    print('Hertz guess:', A, d0)

    i1 = (y < Fnoise).nonzero()[0]
    if len(i1) > 0:
        i1 = i1[0]
        dx = x[i0:i1]
        dy = y[i0:i1]
        indx = dy > 0
    else:
        indx= i1

    if indx.sum() < 4:
        print('Unable to test exponential fit, opening range %d to %d' %(i0, len(x)))
        dx = x[i0:]
        dy = y[i0:]
        indx = dy > 0

    if indx.sum() < 4:
        print('Worst case scenario for exp')
        ft = [0, log(abs(y[i0]))]
    else:
        ft = polyfit(dx[indx], log(dy[indx]), 1)
    # exp(ft[1]) * exp(ft[0]*x)
    print('Exp guess:', exp(ft[1]), -ft[0])
    if verbose:
        pl.plot(dx, exp(polyval(ft, dx)), '-', color= 'lightgreen', linewidth= 2)

    if with_i0:
        fit = leastsq(err_H_exp, [A, d0, i0, exp(ft[1]), -ft[0]],
                    (x, y),
                    full_output= True)
        A, d0, i0, amp, b = fit[0]
    else:
        fit = leastsq(err_H_exp_i0, [A, d0, exp(ft[1]), -ft[0]],
                    (i0, x, y),
                    full_output= True)
        A, d0, amp, b = fit[0]
    ############################
    # Evaluate the results:
    y_fitted= f_H_exp(x, A, d0, i0, amp, b)
    l = 1/b if b != 0 else 0

    if verbose:
        pl.plot(x, y_fitted, '-', color= 'brown', linewidth= 2, alpha= 0.5)
        pl.legend(['data','Hertz fit', 'Hertz error', 'noise threshold', 'exp fit', 'full fit'])
        pl.plot(x[i0], y[i0], 'x', color= 'blue')

    infodict = fit[2]
    #infodict['fvec'] is the last output of err_Hertz_cone
    chi2 = (infodict['fvec']**2).sum()
    #reduced residuals:
    chi2red= chi2/(len(x) - len(fit[0])) #len(fit[0]) is the number of parameters

    msg = fit[3]
    #a definition of r squared regression coeff.
    r2 = 1.0 - (chi2 / ((y-y.mean())**2).sum())

    const_H = 3*(1-nu*nu)/(4.0*sqrt(R))
    E = A*const_H

    # do we have a covariance matrix?
    if fit[1] is not None:
        cov = fit[1]*chi2red
        if with_i0:
            dA, dd0, di0, damp, db = sqrt(diag(cov))
        else:
            dA, dd0, damp, db = sqrt(diag(cov))
            di = 0

        dl = db/b**2 if b != 0 else 0
        dE = dA*const_H
        singular = False
    else:
        cov = None
        if with_i0:
            dA, dd0, di0, damp, db = (-1, -1, -1, -1, -1)
        else:
            dA, dd0, damp, db = (-1, -1, -1, -1)
            di = 0

        dl = -1
        dE = -1
        singular = True
    # end filling up errors


    return {'x': x, 'y': y, 'fitted': y_fitted,
            'A': A, 'dA': dA, 'd0': d0, 'dd0': dd0,
            'E': E, 'dE': dE, 'amp': amp, 'b': b,
            'damp': damp, 'db': db, 'l': l, 'dl': dl,
            'Message': msg, 'r2': r2, 'chi2': chi2,
            'cov': cov,
            'singular': singular}
# end fit_H_exp
