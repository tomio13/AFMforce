""" Fir a general JKR model with surface tension onto indentation data
"""
from numpy import abs, sqrt, exp, log, ones, polyfit, polyval, zeros, pi
from scipy.optimize import leastsq
from matplotlib import pyplot as pl
#from  AFMforce import lm
from AFMforce.LM import lm

def f_gamma(x, x0, A, B):
    """ Calculate adhesion force vs. indentation as:
        for all x > x0:
        A (x-x0)**(3/2) + B*(x-x0)
        0 otherwise
        returns the array of force values
    """
    dx = x - x0
    dx = dx*(dx > 0)
    return A*dx**1.5 + B*dx
#end of f_gamma

def J_no_gamma( parms, *args):
    """ Calculate the error array for leastsq
        parms:  x0, A
        args:   x, y, B, (w)
        In this version B is fixed by the user, thus passed
        at the end of x,y, and an optional w

        return:
        an array of [dy/dx0, dy/dA]
    """
    x0, A = parms

    Na = len(args)
    x = args[0]
    B = args[2]

    dx = x - x0

    boundary_weight = -1E5
    boundary_A = boundary_weight if A <= 0 else 0

    dx = dx*(dx > 0)

    dyda = dx**1.5 + boundary_A
    dydx0 = -1.5*A*sqrt(dx) - B
    return [dydx0, dyda ]
#end of J_no_gamma

def J_gamma( parms, *args):
    """ Calculate the error array for leastsq
        parms:   x0, A, B
        args:   x, y, (w)

        return:
        an array of [dy/dx0, dy/dA, dy/dB]
    """
    x0, A, B = parms
    x = args[0]
    dx = x - x0

    boundary_weight = -1E5
    boundary_A = boundary_weight if A <= 0 else 0
    boundary_B = boundary_weight if B <= 0 else 0

    dx = dx*(dx > 0)

    dyda = dx**1.5 + boundary_A
    dydx0 = -1.5*A*sqrt(dx) - B
    dydb = dx + boundary_B
    return [dydx0, dyda, dydb]
#end of J_gamma

def err_no_gamma( parms, *args):
    """ error function for leastsq
        Again, the B is now set by the user as fixed

        Use A >0 and B>0 forced by a constain employing
        a boundary weight multiplier to penalize the error.
        This still may accept a small negative value,
        lower than -1/boundary_weight.

        parms are x0, A
        arguments: x, y, B and w if it is specified
        returns:
        an array of yfit - y or w*(yfit-y) values
    """
    x0, A = parms
    x = args[0]
    y = args[1]
    B = args[2]

    if len(args) > 3:
        w = args[3]
        noweight = False
    else:
        w = ones( x.shape, dtype='float')
        noweight = True
    #end if there are weights

    boundary_weight = 1E5
    #max(-A,0) will be 0 if A > 0, else something
    #max(-B,0) is similar. These will pull the error up...
    b_A = -A if A < 0 else 0
    boundary_chi = boundary_weight * b_A
    dx = x - x0
    dx = dx*(dx > 0)

    yy = A*dx**1.5 + B*dx
    chi = (yy-y)  if noweight else (yy-y)*w
    return chi+boundary_chi
#end of err_no_gamma

def err_gamma( parms, *args):
    """ error function for leastsq
        parms are x0, A, B
        arguments: x, y, and w if a third one is
        specified
        returns:
        an array of yfit - y or w*(yfit-y) values
    """
    x0, A, B = parms
    x = args[0]
    y = args[1]

    if len(args) > 2:
        w = args[2]
        noweight = False
    else:
        w = ones( x.shape, dtype='float')
        noweight = True
    #end if there are weights

    boundary_weight = 1E5
    #max(-A,0) will be 0 if A > 0, else something
    #max(-B,0) is similar. These will pull the error up...
    b_A = -A if A < 0 else 0
    b_B = -B if B < 0 else 0
    boundary_chi = boundary_weight*(b_A + b_B)
    dx = x - x0
    dx = dx*(dx > 0)

    yy = A*dx**1.5 + B*dx
    chi = (yy-y)  if noweight else (yy-y)*w
    return chi+boundary_chi
#end of err_gamma


def fit_gamma(x,y, Fmax= None, R= 1.0, nu= 0.3,\
        gamma = -1,\
        withweight = False, verbose = True):
    """ fit a general a*(x-d0)**1.5 + B*(x-d0) type curve
        on a background corrected data set.

        If withweight set, use a wegith on the y >0 5*sigma values
        Parameters:
            x,y input arrays
            Fmax: if specified, limit to fit y < Fmax
            R:      radius of ball to calculate E modulus
            nu:     Poisson ratio to calculate E modulus
            gamma:  if >= 0 then use it as a fixed gamma, else fit it
            withweight: if set, use weights on points above 5 sigma background
            verbose:    provide  aplot

        return
        a dict containing the fir results
        A,B, d0                 fit paramters
        dA, dB, dd0:            estimated errors from the covariance matrix
        E, Gamma, dE, dGamma    derived parameters
        fit:                    the whole fit tuple from scipy.optimize.leastsq
        Message:                fit message
        x, y, fitted:           fitted x,y data and fit results (y values)
        chi2, chiered and r2:   fit quality
    """
    #we are fitting B = 2 pi gamma, not gamma dirrectly
    Bgamma = 2*pi*gamma
    N = len(x)
    if len(y) != N:
        raise ValueError('Not matching x,y data set')
    #end length check
    indx = ones( x.shape, dtype= 'bool') if Fmax is None else y < Fmax

    #we use the data up to Fmax only:
    xx = x[indx]
    yy = y[indx]

    if N  < 10:
        raise ValueError('Too short data set')
    #end length of valid segment

    if not (xx > 0).any():
        print("All negative x values")
        #putting x0 to the middle
        xx = xx - xx.mean()
    #end if no positive x

    #Problems to solve:
    #fit depends on first guess a lot, we need some good ones
    #and indentation may be way off, thus we need a  good x0
    N = len(xx)
    Ntail = max( 10 , N/10) if N > 100 else N/ 10
    Ntail = int(Ntail)
    #where are valid y values? Estimate tail noise:
    yy_min = yy[:Ntail].mean() if yy[0] < yy[-1] else yy[-Ntail:].mean()
    yy_sd = yy[:Ntail].std() if yy[0] < yy[-1] else yy[-Ntail:].std()
    print("STD:", yy_sd)

    #estimate a polynomial fit on the data:
    plindx = (yy > 5*yy_sd) & (xx > 0)
    plx = xx[plindx]
    ply = yy[plindx]
    x0= plx.min()
    #plfit = polyfit(plx, yy[plindx], 3, rcond=1E-3)
    #plfit = lm( yy[plindx],[plx**1.5, plx]) if gamma < 0 else lm( yy[plindx]-gamma*plx,[plx**1.5])
    plx = plx - x0
    if gamma < 0:
        plfit = lm( ply,[plx**1.5, plx])
    else:
        ply = plx[plx >0]
        plx = plx[ plx > 0]
        #set plB  as well
        plfit = {'b':[((ply - Bgamma*plx) / plx**1.5).mean(), Bgamma] }
    #end if gamma -- get an estimate of A and B for a start

    print("Polynomial fit gave:", plfit['b'])
    print('X0:', x0)

    #A = abs(plfit[0])
    #B = abs(plfit[1])

    plA = plfit['b'][0] if plfit['b'][0] >= 0 else 0
    plB = plfit['b'][1] if plfit['b'][1] >= 0 else 0
    #x0 = 0.0
    if withweight :
        w = ones(yy.shape)
        w[ yy > 5*yy_sd ] = 5

    if gamma >= 0:
        ferr = err_no_gamma
        fJ = J_no_gamma
        fargs = (xx, yy, Bgamma, w) if withweight else (xx, yy, Bgamma)
        fparms = (x0, plA)
    else:
        ferr = err_gamma
        fJ = J_gamma
        fargs = (xx, yy, w) if withweight else (xx, yy)
        fparms = (x0, plA, plB)
    #end setting up fit functions
    print("Start parameters", fparms)

    #fit = leastsq( err_gamma, [x0, plA, plB], (xx, yy, w),\
    fit = leastsq( ferr, fparms, fargs,\
            Dfun= fJ,\
            ftol = 1E-9,\
            col_deriv = 1, full_output = True)

    #let us work on the resutls:
    if gamma >= 0:
        x0, A = fit[0]
        #we do this, so all goes automatically correct
        B = Bgamma
    else:
        x0, A, B = fit[0]

    yf = f_gamma(xx, x0, A, B)

    if verbose:
        print(fit[0])
        print('Number of calls:', fit[2]['nfev'])
        print('MSG:',fit[3])

        pl.clf()
        pl.plot(x,y, 'bo')
        pl.plot(plx, plA*plx**1.5 + plB*plx, 'g+-')
        pl.plot(xx, yf, 'r-')

    #error estimation:
    infodict = fit[2]
    chi2 = (infodict['fvec']**2).sum()
    #reduced error: use the number of parameters from the result fit[0]
    chi2red = chi2/(len(xx) - len(fit[0]))
    r2 = 1.0 - chi2 / ((yy-yy.mean())**2).sum()

    if fit[1] is not None:
        singular = False
        cov = fit[1]*chi2red
        dx0 = sqrt(cov[0,0])
        dA = sqrt(cov[1,1])
        dB = sqrt(cov[2,2]) if gamma < 0 else 0
    else:
        singular = True
        dA = dB = dx0 = -1.0

    const = 2.0*(1-nu*nu)/(4.0*sqrt(R))
    E = const*A
    dE = const*dA if not singular else -1.0
    Gamma = B/(2*pi) #mN/m if F is nN, x if micrometer
    dGamma = dB/(2*pi)

    result = {'E': E, 'dE':dE, 'Gamma':Gamma, 'dGamma':dGamma,\
            'A': A, 'dA': dA, 'B': B, 'dB': dB, 'd0': x0, 'dd0': dx0, 'fit': fit,\
            'singular': singular,\
            'r2': r2, 'chi2': chi2, 'chi2red': chi2red, 'Message': fit[3],\
            'x': xx, 'y': yy, 'fitted': yf}
    return result
#end of fit_gamma
