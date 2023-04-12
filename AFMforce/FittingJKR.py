#!/usr/bin/env python
""" Try JKR theory for AFM. Develop the fit algorithm
    This file is the development source for testing purposes.
    The key formulation related \delta^1.5 to the f(F), where
    the problem is that if there is a shift in \delta, we have the
    form:
    (\delta - \delta_0)^1.5 = f(F) = \frac{3 R}{4 E*} (F + F_0 \sqrt{2F_0 F + F_0^2} )

    F_0 = -2*F_{min}
"""

from AFMforce import *
from numpy import zeros, ones, sqrt, pi, diag, linspace, polyfit, polyval
from scipy.optimize import leastsq

def f_JKR(F, A, B, d_0):
    """ calcualte the JKR theory in its simplified form, inverse from force to
        indentation.
        Consider: \delta(F) = (3R/(4E*)(F + 2F_0 + \sqrt(4F_0 F + (2F_0)**2)))**(2/3) + \delta_0
        F_0 = -3/2 \pi \gamma R
        in a simplified form:
        \delta(F) = (A (F + B + \sqrt(2BF + B**2)))**(2/3) + \delta_0
        with 3 parameters

        Parameters:
        F       the forces
        A       or 3 R /(4 E*), 3 R (1-\nu**2)/(4E)
        B       or -2F_0 = 3 \gamma \pi R
        d_0     zero indentation (\delta_0_

        return
        indentation values in micrometer
    """
    # we need to limit all to the contact part
    # above -F_C and above d_0

    # what about the units?
    # F is in nN, then F_C is in nN
    # K is in kPa
    # R is in micron
    # 'a' is automatically in micrometers
    delta = zeros(F.size)
    if B < 0 or A < 0:
        return delta

    F2 = zeros(F.shape)
    indx = F > -B/2
    F2[indx] = B + sqrt(2*F[indx]*B + B**2)
    delta[indx] = (A*(F[indx] + F2[indx]))**(2/3.0) + d_0

    return delta
# end of f_JKR

def err_JKR(parms, *args):
    """ error of fit function needed for leastsq to fit the data.

        Parameters:
        parms:  A, B, d_0
                the role of these see above for f_JKR

        *args   F and delta

        return
        the error between delta and f_JKR
    """
    NA = len(args)
    N = len(parms)

    if NA == 2:
        F = args[0]
        d = args[1]
    else:
        raise ValueError("Invalid arguments for fit!")

    if N == 3:
        A, B, d_0 = parms
    else:
        raise ValueError("Invalid parameters!")

    factor = 0
    if A < 0:
        factor += -A*1E6
    if B < 0:
        factor += -B*1E6

    df = f_JKR(F, A, B, d_0)
    err = (df - d) + factor
    # print(err.sum())

    return err
#end err_JKR

def fit_JKR(delta, F, R= 2.0, nu= 0.3, Fmax = None, N= 200):
    """ Force is in nN, R and delta is in micrometers.
        E is in kPa.
        E* = E/(1-\nu**2), reduced Young's modulus

        Use the form of the general equation as:
        (\delta - \delta_0)**1.5 = 3R(1-\nu)**2/(4E) (F + 2F_0 +
        \sqrt(4*F*F_0 + (2F_0)**2))

        Simplified: A = 3 R (1-\nu**2)/(4E), B = 2F_0
        and the minimum force is -F_0.

        Because this problem can be solved only calculating delta for
        every F values, we can easily see that the zero force part of the
        curve is not usable.
        The segment from the maximal adhesion force is selected up to the
        contact (increasing indentation, delta). The rest is discarded.
        The error depends on the number of points used for the fit, so
        one has to be careful with this selection.

        instead of the original fit data, also return an interpolated curve
        to remove the noise inherit in the force values.

        Based on:
        U.D. Schwarz, J. Coll. Int. Sci. 261: 99-106 (2003)

        Parameters:
        delta   all indentation values
        F       corresponding force
        R       radius of indenter ball
        nu      Poisson ration, 0.3 is good for hydrogels
        Fmax    maximum up to this force
        N       number of interpolated fitted points to return

        Return
        a dict containing:
        E*           related to the Young's modulus
        E
        nu
        gamma       from B = 3 \gamma \pi R
        d_0         surface contact point
        F_C         maximal adhesion force (global minimum of the fitted F)
        F, delta:   the fitted segments (original data)
        fitted      the fitted result
        Fi, di      interpolated fitted results
    """
    F_C= - F.min()
    # first threshold is a maximal force:
    indx = F < Fmax if Fmax is not None else ones(F.shape, dtype='bool')

    if indx.sum() < 3:
        raise ValueError('Invalid force limit %.3F' %Fmax)

    xx = F[indx]
    yy = delta[indx]
    # pick where the curve is minimum, we have our first guess
    # for \delta_0
    d_0 = yy[(xx == xx.min())][0]
    B = 2*F_C
    # which part of the curve to use?
    # if x is indentation, then:
    print(d_0, F_C)
    dindx = yy > d_0

    if dindx.sum() < 3:
        raise ValueError('No indentation is left')
    xx= xx[dindx]
    yy= yy[dindx]

    # guess A:
    xa = xx + B + sqrt(2*B*xx + B**2)
    ft = polyfit(xa, (yy - d_0)**1.5, 1)
    #pl.clf()
    #pl.plot(xa, (yy-d_0)**1.5)
    #pl.plot(xa, polyval(ft, xa), 'r-')
    A = ft[0] if ft[0] > 0 else 1.0

    #print('start with:', A, B, d_0)
    # use without explicit Jacobian:
    fit = leastsq(err_JKR, [A, B, d_0], (xx, yy),\
            # Dfun = J_JKR,
            ftol= 1E-8,
            col_deriv= 1,
            full_output= True)

    A, B, d_0= fit[0]
    infodict= fit[2]

    #print('Fit:', fit[0])
    #print('info:', fit[3])

    # the last error output:
    chi2 = (infodict['fvec']**2).sum()

    # the reduced error is divided by the degree of freedom:
    chi2red = chi2/(len(xx) - len(fit[0]))
    msg = fit[3]

    # correlation coefficient squared:
    r2 = 1.0 - (chi2/((yy-yy.mean())**2).sum())

    # do we have a covariance matrix?
    # then estimate the error (inferior)

    if fit[1] is not None:
        dA, dB, dd_0 = sqrt(diag(fit[1]*chi2red))
        singular= False
    else:
        dA, dB, dd_0= (-1, -1, -1)
        singular= True
    #end if

    fitted = f_JKR(xx, A, B, d_0)
    # number of meaningfully fitted points:
    Nf = (xx > -B/2).sum()

    Est = 3*R/(4*A)
    dEst = 3*R*dA/(4*A**2) if dA > 0 else -1
    E = Est*(1-nu**2)
    dE = dEst*(1 - nu**2) if dEst > 0 else -1
    gamma = B/(3*pi*R)
    dgamma = dB/(3*pi*R) if dB > 0 else -1
    # B is positive
    F_C = B/2
    dF_C = dB/2 if dB > 0 else -1
    # interpolated curve for a clear plot of the trend:
    xi = linspace(xx.min(), xx.max(), N)
    yi = f_JKR(xi, A, B, d_0)

    res = {'F': xx, 'delta': yy, 'fitted': fitted,
           'FC': F_C, 'E*': Est, 'E': E, 'gamma': gamma, 'd0': d_0,
           'dFC': dF_C, 'dE*': dEst, 'dE': dE, 'dgamma': dgamma, 'dd0': dd_0,
           'chi2': chi2, 'r2': r2, 'Nfitted': Nf,
           'Fint': xi, 'dint': yi,
           'message': msg}
    return res
#end fit_JKR
