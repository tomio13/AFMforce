#!/usr/bin/env python
""" fit a sine wave background to AFM data and try guessing the parameters
"""

from AFMforce.Smooth import Smooth
from matplotlib import pyplot as pl
from numpy import sin, pi, diff, polyfit, polyval, abs, log, floor
from scipy.optimize import leastsq

__all__ = ['bg_wave_fit', 'fit_wave_background']

def bg_wave(x,
                amplitude= 1.0,
                amplitude_slope= 0.0,
                omega= 0.1,
                beta= 0.0,
                delta= 0.0,
                offset= 0.0,
                slope= 0.01):
    """ a wavy background function to handle interference background
    """
    y = (amplitude + amplitude_slope*x)* sin((omega + beta*x) * x + delta) + offset + slope*x

    return(y)
# end bg_wave


def err_bg_wave_fit(params, *args):
    x = args[0]
    y = args[1]
    amplitude, amplitude_slope, omega, beta, delta, offset, slope = params
    # amplitude = abs(amplitude)
    amplitude = amplitude

    return (y - bg_wave(x, amplitude, amplitude_slope, omega, beta, delta, offset, slope))
# end calculating error vector


def fit_wave_background(distance, force,
                        smooth_radius= 30,
                        verbose= False):
    """ fit a sine wave background using a 7 parameter formula,
        accounting for: constant and linear offset,
        phase, and chirp,
        linear amplitude shift.

        The algorithm first subtracts a linear fit from the data,
        thus make sure you have enough oscillations in them to get a
        reasonable cut through.

        Then it determines the spatial frequency analyzing the zero
        crossing of the linear corrected data, smoothed with Smooth
        and using the smooth_radius.

        Parameters:
        distance:       distance or Z position
        force:          deflection of the cantilever
        smooth_radius:  smooth width passed to Smooth
        verbose:        plot the fit result


        Return a dict with:
        {'distance': the x-array provided
         'fitted': the fitted data
         'fit': the full response from leastsq
         }
    """
    if len(distance) != len(force):
        raise ValueError('length mismatch in arrays!')

    # guess order of X-axis:
    order_mag = 10**(floor(log(distance.max() - distance.min())/log(10.0) + 0.5))

    # get a linear fit to the background as a guess
    fit_line = polyfit(distance, force, 1)
    force_corr = force - polyval(fit_line, distance)

    # positive to negative transition in force happens at every half
    # wavelength... to guess omega, we use that:
    # omega_guess = pi/diff(distance[diff(Smooth(force_corr, 30) > 0).nonzero()[0]]).mean()
    # try underestimating it:
    omega_guess = pi/diff(
            distance[diff(Smooth(force_corr, smooth_radius) > 0).nonzero()[0]]
            ).mean() /2.0
    # print('guessing omega:', omega_guess)
    # print('amplitude guess:', force_corr.max())


    # now, get a background fit
    bg = leastsq(err_bg_wave_fit,

                 [force_corr.max(), -0.001/order_mag,
                  omega_guess, 0.001/order_mag, 0.0,
                  fit_line[1], fit_line[0]],

                 (distance, force),
                 maxfev= 1000,
                 col_deriv= True,
                 full_output= True)

    fitted_force = bg_wave(distance,
                         bg[0][0],
                         bg[0][1],
                         bg[0][2],
                         bg[0][3],
                         bg[0][4],
                         bg[0][5],
                         bg[0][6])

    if verbose:
        pl.plot(distance, force, alpha= 0.5)
        pl.plot(distance,
             fitted_force,
             'r-',
             linewidth= 2)
        pl.xlabel('Z position')
        pl.ylabel('vertical deflection')

    res =  {'distance': distance,
            'force': force,
            'fitted': fitted_force,
            'fit': bg
            }
    # to play in the fit parameters with their names, we take the names
    # and the values from bg[0]:
    keys = ['amplitude', 'amplitude_slope', 'omega', 'beta', 'delta', 'offset', 'slope']

    res.update({i:j for i,j in zip(keys, bg[0])})

    return(res)
# end fit_wave_background
