#!/usr/bin/env python
""" Identify and fit break off forces allowing for classification
    of the events.
    Here we collect funcitons to fit the various rip-off types using:
    worm-like chain model
"""
from numpy import abs, sqrt, polyfit, polyval, isnan, zeros, ones
from numpy.linalg import pinv
from scipy.optimize import leastsq
from matplotlib import pyplot as pl
kB = 1.3806504 #the Boltzmann constant, without its 1E-23 factor

###############################################################################
# WLC fitting

def f_WLC(x, T, lc, L, F0= 0.0, Lunit=1E-6, Funit=1E-9):
    """ Calculate the force according to the worm-like chain model
        F = kB T/lc (1/(4*(1-x/L)^2) + x/L - 1/4)
        This force is only valid in the 0 <= x < L interval
        We set it to 0 otherwise (should we set nan...)

        Parameters:
            x       position data set
            T       temperature
            lc      persistence length in micrometers
            L       contour lengh
            F0      constant force background
            Lunit   unit of length, 1E-6 micron
            Funit   unit of force, 1E-9 nN

        return:
            F       an array of forces to each x value in nN
    """
    if lc <= 0:
        raise ValueError("Invalid persistence length (<= 0)!")

    #assuming L and x have the same dimensions, the sequence in the
    #parenthesis is dimensionless.
    #lc should be in microns = 1E-6 m
    #F then in nN, 1E-9 N
    #thus 1E-23 / 1E-6 = 1E-17; 1E9 * 1E-17 = 1E-8
    kBfactor = 1E-23/(Lunit*Funit)

    #anyone throwing invalid parameters at us?
    if isnan(lc) or lc <= 0 or isnan(L) or L==0:
        return zeros(x.shape)

    xx = x
    F= kB*kBfactor*T/lc *( xx/L - 0.25 + 1.0/(4.0*(1-xx/L)**2)) + F0
    F = F*(xx < L)
    return F
#end of f_WLC

def J_WLC(parms, *args):
    """ Estimate the values of the Jacobian.
        Based on C*(x/L - 0.25 + 1/(4*(1-x/L)**2)

        Returns:
        [dF/dC, dF/dL]
    """

    x = args[0]
    N = len(parms)
    if N == 3:
        C, L, F0 = parms

    elif N == 2:
        C, L = parms
        F0 = 0.0
    else:
        raise ValueError("Invalid set of parameters!")

    if isnan(L) or L==0:
        #print("Invalid parameters!")
        #for invalid parameters the whole landscape is 0...
        #the retusn array is either 2 or 3, shape of x:
        shape = (N,) + x.shape
        return zeros(shape)

    #our landcape is only valid for 0 < x < L
    #if x > L -> we have overstretched the molecule, the model is gone
    indx = x < L

    L2 = L*L
    dC = zeros(x.shape)
    dL = zeros(x.shape)
    dx0 = zeros(x.shape)
    #shortened version of dx values:
    dxs = x[indx]
    dxL3 = (1 - dxs/L)**3

    dC[indx] = dxs/L - 0.25 + 1.0/(4.0*(1-dxs/L)**2)
    dx0[indx] = -1.0*C/L*(1.0 + 1.0/(2*dxL3))
    #dL[indx] = -1.0*C*dxs/L2*(1.0 + 1.0/(2*dxL3))
    dL[indx] = dxs/L * dx0[indx]
    dF0 = ones(x.shape)

    #print('J_WLC with',N, 'parameters', len(dL), len(dC))
    if N > 2:
        return [dC, dL, dF0]
    #else

    return [dC, dL]
#end of J_WLC

def err_WLC(parms, *args):
    """ chi^2 function for the WLC model.
        To fit the WLC model, we have a couple of unknowns:
        1. F0 the background force -> we may force it to 0
        2. Because T and lc are coupled, there is only a constant
            prefactor in the form:
            C*(x/L - 0.25 + 1/(4.0*(1-x/L)**2)
            for x > 0, otherwise 0


        parms: C, L
        args: x, y, weight

        Return:
        (yy-y)*weight or (yy-y) if weight is not provided

        (leastsq does the square and sum up)
    """
    NA = len(args)
    N = len(parms)

    if NA > 3 or NA < 2:
        raise ValueError("Invalid Parameters")

    if N == 3:
        C, L, F0 = parms
    elif N == 2:
        C, L = parms
        F0 = 0.0
    else:
        raise ValueError("Invalid number of parameters!")

    x = args[0]
    y = args[1]

    indx = (x >= 0)*(x < L)

    if isnan(L) or L <= 0 or C <= 0 :
        #print("invalid guess...")
        return abs(y) * 1000000000

    yy = zeros(x.shape)
    yy[x < 0] = F0
    dxs = x[indx]
    yy[indx] = C*(dxs/L - 0.25 + 1.0/(4.0*(1.0 - dxs/L)**2))+F0

    if NA == 3:
        weight = args[2]
        noweight = False
    else:
        noweight = True

    chi = yy-y if noweight else (yy-y)*weight
    chi[x > L] = 100 * y.max()

    #print (chi**2).sum()
    return chi
#end err_WLC

def fit_WLC(x,y, T=298.0, pars=[], \
        Lunit = 1E-6, Funit=1E-9,\
        verbose= False):
    """ fit the WLC model to a data set.

        Parameters:
        x,y     force distance curve, possibly x in micron, y in nN
        T       temperature of the measurement (for lc)
        pars    optional start parameters
                if provided, it has to be [Fmax, L]
        Lunit   the unit of L: 1E-6 micron, 1E-9 nm
        Funit   the unit of force: 1E-9 nN, 1E-12 pN

        return:
        a dict containing
        "z", "F":   the original x,y data
        "Ffitted"   fitted force values
        "lc"        persistence length
        "T"         temperature used
        "L"         countour length
        'Ft'        thermal driving force, kB T/lc
        "fit"       all fit parameters from leastsq
        ...
    """
    N = len(x)
    kBfactor = 1E-23 / (Lunit*Funit)

    if N > 3 and N != len(y):
        raise ValueError("We need > 5 data points and equal arrays!")

    #end if length check

    #try a simple estimate if the force curve goes to negative or
    #positive... If negative, turn it around!
    if abs(y.min()) > abs(y.max()):
        if verbose:
            print("negating y")
        y = -1.0*y

    if len(pars) == 2:
        C, L = pars

    else:
        #rough estimate of the fitting parameters:

        #the force is calculated as:
        #    F= kB*1E-8*T/lc *( xx/L - 0.25 + 1.0/(4.0*(1-xx/L)**2))
        #the prefactor is in the order of:
        C = 0.01*(y.max() - y.min())
        L = 1.05*x.max()

    if verbose:
        print("Start values C: %.5f, L: %5f" %(C, L))

    #fit the rest:
    fit = leastsq( err_WLC, [C, L], (x, y),\
            Dfun = J_WLC, \
            #relative tolerance for the sum of squares:
            ftol = 1.0E-12,\
            #relative tolerance desired for the solution
            xtol = 1.0E-12,\
            #maximal number of iterations:
            maxfev=10000,\
            col_deriv=True,\
            full_output = True)

    #####################################################################
    ######### RESULTS:
    C, L = fit[0]
    info = fit[2]
    if verbose:
            print("End values C: %.5f, L: %.5f" %(C,L))

    #result:
    #calculate the lc back from the constant:
    if C == 0.0:
        print("Null force!")
        lc = 0.0
    else:
        lc = kB*T*kBfactor/C

    #f_WLC handles F0 as an additive constant, no problem:
    yy = f_WLC(x, T, lc, L, 0, Lunit, Funit) if lc != 0.0 else zeros(y.shape)

    #info['fvec'] is the last output of the model
    chi2 = ((yy-y)**2).sum()

    #reduced residuals:
    #division with the number of freedom
    chi2red = chi2/(len(x) - 2) #2 is the number of parameters to the fit
    #r-square, the regression coefficient 1 - RSS/TSS
    r2 = 1.0 - (chi2/((y-y.mean())**2).sum())

    if fit[1] is not None:
        cov = fit[1] * chi2red
        dC = sqrt(abs(cov[0,0]))
        dL = sqrt(abs(cov[1,1]))
    else:
        print("singular value in fit!")
        cov = None
        dC = -1.0
        dL = -1.0
    #end if Hessian defined

    dlc = lc * dC/C if C != 0 else 0.0

    if verbose:
        pl.clf()
        pl.plot(x,y,'go', alpha=0.3)
        pl.plot(x, yy, 'r-', alpha=0.3)

#    return {'Ft': C, 'L': L, 'lc': lc, 'F0': F0,\
#            'dFt':dC, 'dL':dL, 'dF0':dF0, 'dlc':dlc, 'T': T,\
#            'fit': fit, 'Z':x, 'F':y, 'fitted': yy, 'r2':r2, 'chi2':chi2,\
#            'cov':cov, 'Message': fit[3]}
    return {'Ft': C, 'L': L, 'lc': lc,\
            'dFt':dC, 'dL':dL, 'dlc':dlc, 'T': T,\
            'fit': fit, 'Z':x, 'F':y, 'fitted': yy, 'r2':r2, 'chi2':chi2,\
            'cov':cov, 'Message': fit[3]}
#end of f_WLC
