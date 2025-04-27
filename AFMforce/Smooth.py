#/usr/bin/env/python
""" Smoothing functions for AFM curves
"""
import numpy as nu
from numpy.fft import fft, ifft

######## Smoothening:
def SavGol(y, Nkernel, order=4, difforder=0, KernelOnly=False):
    """ Implement a Savitsky-Golay smoothing of a 1D data vector.
        Based on the savgol filter in R by Hans W. Borchers (2003),
        and the example in the nympy cookbook by Ryan Hope
        (https://gist.github.com/RyanHope/2321077).
        Due to issues with the padding / shifting in nu.convolve,
        implement 1D FFT based convolution.

        Further references:
            A. Savitzky, M. J. E. Golay, Smoothing and Differentiation of
            Data by Simplified Least Squares Procedures. Analytical
            Chemistry, 1964, 36 (8), pp 1627-1639.

            Numerical Recipes in C++
            3rd Edition
            Press, William H.; Teukolsky, Saul A.;
            Vetterling, William T.; Flannery Brian P.

            Cambridge University Press ISBN-13: 9780521880688

        Parameters:
        y           dataset
        Nkernel     width of the kernel
        order       order of smoothing
        difforder   the order of differential to be returned
    """
    N = len(y)
    Nkernel = int(Nkernel)
    if Nkernel < 1:
        raise ValueError("Too small kernel window!")

    if order < 1:
        raise \
        ValueError("Too small order, the fit should be at least linear (1)")
    #we go half the kernel +/- the midpoint
    #integer will truncate it:
    Nk2 = int(Nkernel/2)
    #we need an outer product using the power of the original data
    #0 to order including both ends!
    #the numpy.outer can not do this, though it were faster
    #M = mat( [[k**i for i in range(order)] for k in range(-Nk2, Nk2+1)])
    #let us use a bit more of the numpy:
    #this is faster for large numbers, but with a small margin only
    # M = nu.mat( [nu.arange(-Nk2,Nk2+1)**i for i in range(order+1)]).T
    M = nu.asmatrix( [nu.arange(-Nk2,Nk2+1)**i for i in range(order+1)]).T
    #solve, and get the line of the order of derivative:
    Mi = nu.linalg.pinv(M).A[difforder]
    #up to this point it matches the results of R
    if KernelOnly:
        return Mi

    #do the convolution using fft:
    lg2 = nu.log(2)
    #round length up to the next power of two
    Nlog2 = 2**int(nu.log(N+Nkernel)/lg2+1)
    Ndiff = int((Nlog2 - N)/2)
    #y = nu.concatenate( (y, nu.zeros(Nlog2 - N)+ y[-1]) )
    #use symmetric padding to remove edge jumps:
    y = nu.concatenate( (nu.zeros(Ndiff) + y[0], y,\
                         nu.zeros(Nlog2 - N - Ndiff)+ y[-1]) )

    res = ifft( fft(y, Nlog2)*fft(Mi[::-1], Nlog2))

    Nstart = Ndiff + Nk2
    #to make sure the shift is a match!
    #return (-1)**difforder*res[Nk2:N+Nk2].real
    return (-1)**difforder*res[Nstart:N+Nstart].real
#end of SavGol()



def Smooth(f, radius=10, par=2, type="SavGol", deriv= 0):
    """ Smoothen the data using a convolution filter of
        a Gauss kernel or a Savitsky-Golay kernel
        Make sure you chose parameters such that the kernel vanishes
        at the edges, or it will shift your curve.
        The curve is padded by its end values to compensate its shift,
        but doing the same with the kernel would kill e.g. a boxcar
        smoothing.

        Parameters:
        radius      radius of the filter. The whole is 2r+1 wide
        par         is sigma for the Gaussian kernel
                    or the order of the fit for the Sav. Gol. filter
                    (2 is a good one for the latter)
        type        what kernel to use.
                    It can be a Gaussian, a boxcar for a flat kernel,
                    or SavGol for the Savitsky-Golay filter (default)
                    or one of the window functions defined in numpy:
                    Bartlett, Hanning, Kaiser or Blackman
        deriv       Optional parameter, calculate the first derivative
                    if it is nonzreo, for the case of Gaussian or
                    Savitzky - Golay filters
                    Note: The Gaussian gives the negative of a real derivative
                    when used on a signal.

        Return:
        the smoothened vector with the same length as f
        Do not forget, that the first and last r points are
        affected by the finite sample length!
    """
    if radius < 2:
        raise ValueError("Window should be wider as 1 point!")
    if deriv != 0:
        deriv = 1

    Nk = int( 2*radius + 1)
    kx = nu.arange(Nk, dtype="f") -radius
    type = type.lower()

    if "gauss" in type:
        if par < 0:
            print("width < 0: using absolute value")
            par = -1.0 * par

        width2 = -2.0*par*par
        #This part is the same:
        k = 1.0/(nu.sqrt(2.0*nu.pi)*par) * nu.exp((kx*kx)/width2)
        if deriv != 0:
            # derivative of the Gaussian kernel is:
            # k = -k*kx/par**2
            # this will cause a negative derivative, the leading part
            # has the positive weight!
            k = k*kx/par**2
        #if it is too narrow, we may have to renormalize it to avoid
        # shifting the data due to offset:
        if radius < 4*par:
            k = k - k.min()
            k = k / k.sum()

    elif "flat" in type or "boxcar" in type:
        k = k / float(k.sum())

    elif "SavGol" in type or "savgol" in type:
        k = SavGol(f, Nk, par, deriv, KernelOnly=True)

    else:
        try:
            txt = "nu.%s(%d)" %(type, Nk)
            k= eval(txt)
            k = k/k.sum()
        except:
            print("Invalid kernel type!")
            return f
    #do the convolution using fft:
    N = len(f)
    Nk2 = int(radius)
    lg2 = nu.log(2)
    #round length up to the next power of two
    #to make sure the shift is a match!
    #we want convolution, so pad the data with zeros, actually
    #lift this to the last point...
    Nlog2 = 2**int(nu.log(N+Nk)/lg2+1)
    Ndiff = int((Nlog2 - N)/2)
    #ff= nu.concatenate( (f, nu.zeros(Nlog2 - N)+ f[-1]) )
    #use symmetric padding to compensate edge jumps
    y = nu.concatenate( (nu.zeros(Ndiff) + f[0], f,\
                         nu.zeros(Nlog2 - N - Ndiff)+ f[-1]) )

    #res = ifft( fft(f, Nlog2)*fft(k[::-1], Nlog2))
    #Subtracting the mean of the signal should improve precision
    # and help removing offsets
    ymean = y.mean()
    res = ifft( fft(y-y.mean(), Nlog2)*fft(k[::-1], Nlog2))
    if deriv == 0:
        #we add the offset back for smoothing:
        res = res + ymean
    #end if offset was done

    #mode "valid": only where the signals overlap
    Nstart = Ndiff + Nk2
    return res[Nstart:N+Nstart].real
#end of Smooth
