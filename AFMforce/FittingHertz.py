#!/usr/bin/env python
""" Fitting functions to be used in relation to AFM force curves.
"""

#from JPKforce import *

from numpy import abs, sqrt, pi, tan, ones, polyfit, polyval
from scipy.optimize import leastsq
from matplotlib import pyplot as pl

def f_Hertz_ball(x, E, R, nu):
    """ Calculate the force for a Hertz model using a spherical
        indenter. Only for x > 0, others are zero.

        Parameters:
            x       indentation data (distance)
            E       Young's modulus
            R       radius of the sphere
            nu      Poisson ratio

        Return
        array of force values for each x

        base on:
        Adv. Coll. Int. Sci. 120: 57-67 (2006)
    """
    #formula: 4/3 sqrt(R) E/(1-nu**2) * x**(3/2)
    F = 4.0*E*sqrt(R)/(3.0*(1.0 - nu*nu)) if nu != 1 else 0
    xn = x*x*x*(x > 0)

    return F*sqrt(xn)
#end f_Hertz_ball

def f_Hertz_cone(x, E, theta, nu):
    """ Calculate the force for a Herz model using a conical indenter
        Only x > 0 has a meaning, all others are zero.

        Parameters:
            x       indentation data
            E       Young's modulus (Pa)
            theta   half opening angle of the cone
            nu      Poisson's ration (most often 0.5 for cells)

        Return:
        array of force values for each x
    """

    F = 2.0*E*tan(theta)/(pi*(1-nu*nu)) if nu != 1.0 else 0.0
    xn = x*(x > 0)
    return F*xn*xn
#end f_Hertz_cone

def J_Hertz_cone(parms, *args):
    """ Calculate the Jacobian for err_Hertz_cone, based on
        A = 2.0*E*tan(theta)/(pi*(1-nu*nu))
        and: F = A*(x-d0)**2 + F0

        parms: x values (indentation)
        args:   A, d0

        Return an array [dF/dA, dF/d d0, dF/dF0] for each x value
    """
    A, d0, F0 = parms
    x = args[0]
    x = x - d0

    dA = x**2
    dd = -2.0*A*x
    dx = (x > 0) if not (x<0).all() else 1000*ones(x.shape, dtype='float')
    #F0 derivative is 1, independent of all others
    return [dA*dx, dd*dx, ones(x.shape,dtype='float') ]
#end of J_Hertz_cone

def J_Hertz_ball(parms, *args):
    """ Calculate the Jacobian for err_Hertz_ball, based on
        A = 4.0*E*sqrt(R)/(3*(1-nu*nu))
        and: F = A*(x-d0)**1.5 + F0

        parms: x values (indentation)
        args:   A, d0

        Return an array [dF/dA, dF/d d0, dF/dF0] for each x value
    """
    A, d0, F0 = parms
    x = args[0]
    x = x - d0

    dx = x*(x > 0) if not (x < 0).all() else 1000*ones(x.shape, dtype='float')
    dA = dx**1.5
    dd = -1.5*A*sqrt(dx)
    #F0 derivative is 1, independent of all others
    return [dA, dd, ones(x.shape,dtype='float') ]
#end of J_Hertz_ball

def err_Hertz_ball(parms, *args):
    """ Error function for fitting f_Hertz_ball.
        parms: A, d0, F0
            where A can be translated as:
            4.0*E*sqrt(R)/(3*(1-nu*nu))
            function: A*x**(1.5) + F0

        data (*args): x, y, weight

        return:
            (fitted - y)*weight

        leastsq will take the square of this and sum it
    """

    NA = len(args)
    N = len(parms)

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

    #if parameters do not match, let it crash:
    A, d0, F0 = parms
    #only positive indentation matters:
    x = x - d0
    x = x*(x > 0)

    #save a call, use the same code as above:
    if A == 0.0 or (x==0).all():
        return 10000*y

    yy = A*x**1.5 + F0

    #but calculate the difference for all points:
    chi = (yy-y)*weight if not noweight else (yy-y)

    return chi
#end err_Hertz_ball

def err_Hertz_cone(parms, *args):
    """ Error function for fitting f_Hertz_cone.
        parms: A, d0, F0
            where A can be translated as:
            2.0*E*tan(theta)/(pi*(1-nu*nu))

        data (*args): x, y, weight

        return:
            (fitted - y)*weight

        leastsq will take the square of this and sum it
    """

    NA = len(args)
    N = len(parms)

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

    #if parameters do not match, let it crash:
    A, d0, F0 = parms
    #only positive indentation matters:
    x = x - d0
    x = x*(x > 0)

    #save a call, use the same code as above:
    if A == 0.0 or (x == 0).all():
        return 10000*y

    yy = A*x*x + F0

    #but calculate the difference for all points:
    chi = (yy-y)*weight if not noweight else (yy-y)

    return chi
#end err_Hertz_cone


def fit_Hertz(x,y, Fmax = None, theta=0.349065885,
            R=2.25, nu=0.5, type="cone", withweight=True, verbose=False):
    """ Fit the Hertz model of a cone to the data x, y,
        where x is indentation, y is the force. Fmax is a
        maximal indentation, below which the values are considered
        all above are ignored.

        Parameters:
        x,y     data set
        Fmax    if specified, use force values below this limit
        theta   the half opening angle of the cone, in radians
        nu      is the Poisson ratio
        R       radius of the partcle (ball) in micrometers
        type    "cone" or "ball" for which type of model to fit
                if "ball", then specify R instead of theta
        withweight  if True, use the abs(y) as weight in fitting
        verbose:    print starting guess and plot a figure
                    points are the raw data, red line is the
                    first guess, green is the fitted result

        return:
        a dict describing the fit

        Note:
        if F was in nanonewtons, x in microns and R in microns,
        then the resulted E has to be multiplied with 1000 to get
        in Pa.
    """

    #filter according to Fmax:
    #if the array is empty, it will cause some trouble, but for now.
    #we consider that user failuer
    #indx = y < Fmax if Fmax != 0 else ones(y.shape, astype='bool')
    indx = y < Fmax if Fmax is not None else ones(y.shape, dtype='bool')
    if indx.sum() < 3:
        raise ValueError("Invalid force limit: %.3f" %Fmax)

    xx = x[indx]
    yy = y[indx]

    #start values:
    if not (xx > 0).any():
        print("all negative or zero X values!")
        xx = xx - xx.mean()
    #else we can go for the positive x values for estimating
    #our guesses:
    xn = xx*(xx>0)

    if type.lower() == "cone":
        typ = 1
    elif type.lower() == "ball":
        typ = 2
    else:
        raise ValueError("Invalid type specification: %s", type)
    #end if

    #get a rough estimate for minima and its noise:
    i_N = len(yy)
    i_min = int(i_N/10)
    yy_min = yy[:i_min].mean() if yy[0] < yy[-1] \
                                    else yy[-i_min:].mean()
    yy_sd = yy[:i_min].std() if yy[0] < yy[-1] \
                                    else yy[-i_min:].std()

    #This is way too much:
    #w = (xx - xx.min())*abs(yy - yy_min)
    #add more weight to those standing out of the 2nd (95%)
    w = (abs(yy-yy_min) > 2*yy_sd) + 1
    #print("Wmin: %.4f, Wmax: %.4f" %(w.min(),w.max()))

    if typ == 1:
        #ft = polyfit(xn, yy*(xx>0),2, w = abs(yy*(xx>0)))
        ft = polyfit(xn, yy*(xx>0),2, w = w)
        A0 = abs(ft[0])
        d0 = abs(ft[1])/(2.0*A0) if A0 != 0 else sqrt(abs(ft[2]))
        F0 = ft[2] - A0*d0*d0

        yfitted = polyval(ft,xn)
        err_f = err_Hertz_cone
        J_f = J_Hertz_cone

    elif typ == 2:
        #make a linearized model using y**2:
        #this should be a 3rd order polynomial containing A**2...
        yfit = yy*yy #only for the polynomial fit!
        #ft = polyfit(xn, yfit*(xx >0), 3, w=abs(yy))
        ft = polyfit(xn, yfit*(xx >0), 3, w= w)
        A0 = sqrt(abs(ft[0]))
        #too many possibilities, anyway this is an approximation only:
        d0 = - ft[1]/(3.0*abs(ft[0])) if ft[0] != 0.0 else abs(ft[2]/ft[1])
        #in the firts approximation: F0 = sqrt(-A**2 * d0**3)
        #F0 = abs( A0*A0*d0)**1.5
        #F0 = ft[3]-ft[2]*d0/3.0
        #It may be better using yy_min
        F0 = yy_min

        yfitted2 = polyval(ft,xn)
        yfitted = sqrt((yfitted2 > 0)*yfitted2)
        err_f = err_Hertz_ball
        J_f = J_Hertz_ball

    #end if typ...
    #for fitting with weight
    #negative or mostly negative force values should be no indentation
    #if weights are used, we kill those weights
    w[yy < (yy_min-2*yy_sd) ] = 0

    #it may be just too close to 0:
    if abs(F0) < 1E-3:
        F0 = 0.0

    if verbose:
        pl.clf()
        pl.plot(xx,yy,'b+',alpha=0.3)
        pl.plot(xn, yfitted,'r-')

        print("polyfit", ft)
        print("Start values (A0,d0,F0)", A0, d0, F0)

    #now do fitting:
    if withweight:
        #fit = leastsq( err_f, [A0,d0, F0], (xx,yy,abs(yy)), \
        fit = leastsq( err_f, [A0,d0, F0], (xx,yy,w), \
            Dfun = J_f, ftol= 1.0E-8,\
            col_deriv =1, full_output=True)
    else :
        fit = leastsq( err_f, [A0,d0, F0], (xx,yy), \
            Dfun = J_f, ftol= 1.0E-8,\
            col_deriv =1, full_output=True)

    A, d0, F0 = fit[0]

    infodict = fit[2]
    #infodict['fvec'] is the last output of err_Hertz_cone
    chi2 = (infodict['fvec']**2).sum()
    #reduced residuals:
    chi2red= chi2/(len(xx) - len(fit[0])) #len(fit[0]) is the number of parameters

    msg = fit[3]
    #a definition of r squared regression coeff.
    r2 = 1.0 - (chi2 / ((yy-yy.mean())**2).sum())

    if typ == 1:
        # A = 2.0*E*tan(theta)/(pi*(1-nu*nu))
        const = pi*(1-nu*nu)/(2.0*tan(theta))
    elif typ == 2:
        # A = 4.0*E*sqrt(R)/(3*(1-nu*nu))
        const = 3.0*(1-nu*nu)/(4.0* sqrt(R))
    else:
        raise ValueError("Invalid type!!!")
    #end if typ
    E = A*const

    #error estimates:
    if fit[1] is not None:
        cov = fit[1]*chi2red
        dA = sqrt(cov[0,0])
        dd0 = sqrt(cov[1,1])
        dF0 = sqrt(cov[2,2])
        singular = False

    else:
        print("Singular value!\n")
        cov = None
        dA = -1
        dd0 = -1
        dF0 = -1
        singular = True
    #end if

    if typ == 1:
        yf = A*(xx > d0)*(xx - d0)**2 + F0
    elif typ == 2:
        ftemp = (xx > d0)*(xx - d0)
        yf = A*ftemp**1.5 + F0

    if verbose:
        pl.plot(xx,yf,'g-')
        pl.draw()

    res = { 'A': A, 'd0': d0, 'F0': F0, 'fit': fit, \
            'dA': dA, 'dE': dA*const, 'dd0':dd0, 'dF0':dF0,\
            'E': E, "singular": singular, "fitted":yf,\
            'Message':msg, 'r2': r2, 'chi2': chi2, "cov":cov,
            'x':xx, 'y':yy}

    return res
#end fit_Hertz_cone

