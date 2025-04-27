#!/usr/bin/env python
""" Scripts to import JPK AFM force data from force spectroscopy files
"""

import os
from zipfile import ZipFile
import numpy as nu
from matplotlib import pyplot as pl
from AFMforce.WavyBgCorr import bg_wave, fit_wave_background
from AFMforce.Smooth import SavGol, Smooth
#from FittingHertz import f_Hertz_cone, f_Hertz_ball, fit_Hertz
#from FittingPull import f_WLC, fit_WLC

#### Routines exported to the external world:
#__all__=["AFMforce", "JPKforce", \
#        "Smooth", "SavGol",\
#        "Baseline", "Contact_point",\
#        "Z_to_D", "find_slope", "FindMinima", "Indentation", \
#        "f_Hertz_cone", "f_Hertz_ball", "fit_Hertz",\
#        "f_WLC", "fit_WLC"]

###############################################################################
########## Actual routines:
class AFMforce(object):
    """
        An abstract class to collect all functions and functional properties of
        AFM force curves. This class is the base for equipment specific
        force classes, such as JPKforce.
        The basic idea is to create an internal dictionary containing:
            header information
            segments
                piexo and or height data, and unit
                time data and unit
                force and deflection data down to sensor level (Volts)
                    and their units
            force constant (if defined)
            sensor response (if defined)
            'segments' -> list of segment names (keys)
            position (X,Y, unit)

        The content of header is optional, we can think about it later
        This class then provides all the information in a uniform manner,
        so analysis can be performed independent of the equipment type.

        When defining a new class for an instument, make sure to modify the
        read_data() function to generate the proper self.dataset and its
        key / value pairs.

    """
    def __init__(self,filename=""):
        """ Initialize the data class and load data if there is anything
            to load.
        """
        self.filename = ""
        self.dataset = {}

        #segment index, could be anything, it is a dict
        #actually used:
        self.__segment__ = -1

        #Z0 is used to set where the surface is for Z -> distance
        self.Z0 = None       #contact point
        #Z.max is the maximal Z value taken into account for analysis
        #Setting Z.max, the data is automatically truncated
        self._Zmax_ = -1.0  #range limit, valid if > 0.0

        #index of data to be returned. Set by setting Zmax
        self._indx_ = nu.empty(0)  #index for filtering

        #perform the actual import:
        if filename != "" and os.path.isfile(filename):
            #perform the import
            self.read_data(filename)

        if self.segments:
            self.segment = self.segments[0]
    #end __init__

    def read_data(self, filename):
        """ The import function of the specific AFM class.
            This function has to take care of filling up the internal
            dict structure with all the keys needed by the class.
            By default this is a dummy function, the specific class should
            implement it!
            All data should be stored in self.dataset
        """
        self.dataset = {}
        #The imported dict has to have a key for each segment
        # e.g. 0, 1, 2 etc.
        # within each segment:
        #   name: informative name, e.g. retract, extend, approach, withdraw, surface
        #
        #   position: x,y, unit -> optional
        #   force: raw_data (in volts), raw_data_unit, deflection (data), deflection_unit
        #       data, unit
        #       if data is not defined, then deflection if not defined, then raw data
        #
        #   height: data, unit, type (measured, piezo)
        #   time: data, unit (time data of the segment)
        #
        #Force constant and sensor response should be global, not segment dependent:
        #   sensor response (in deflection_unit / raw_data_unit)
        #   force constant (in unit / deflection_unit)
        #
        return
    #end import

    @property
    def segment(self):
        return self.__segment__

    @segment.setter
    def segment(self, value):
        """ Set the segment index and reset Z0 as well as check for direction
            The default way of accessing the data segments
        """
        if value in self.segments:
            self.__segment__ = value
            self.Z0 = None

            #invert needs to remove and reapply Zmax
            self.invert()

        else:
            raise ValueError("Invalid segment request %d" %value)
    #end segment

    @property
    def segment_name(self):
        """ The name defined in the headers for this segment"""
        a = self.full_segment

        #does a['name'] exist? From python 3:
        if 'name' in  a:
            return a['name']
        else:
            return ""
    #end segment_name

    @property
    def segments(self):
        """ indices available to be used as segments """
        if "segments" in self.dataset:
            return self.dataset['segments']
        else:
            raise ValueError("No segments defined in dataset!")
    #end def segments

    @property
    def full_segment(self):
        """ a shortcut to the current segment """
        return self.dataset[self.segment]
    #end full_segment

    @property
    def position(self):
        """ Return the [X,Y] position if they are defined.
            Read only parameter
        """
        #if segment is not in a, then it is anyway crappy...
        a = self.dataset

        if "position" in a and "data" in a["position"]:
            return a['position']["data"]
        else:
            return nu.empty(0)
    #end position

    @property
    def position_unit(self):
        """ the unit of the position if defined """
        a = self.dataset

        if "position" in a:
            a = a["position"]
        if "unit" in a:
            return a['unit']
        else:
            return ""
    #end position_unit


    @property
    def force(self):
        """ get a froce curve if it exists
            You can not set this data, but obtain it from
            the raw vertical deflection by setting
            sensor_response and force_constant
        """
        a = self.full_segment

        if "force" in a:
            a = a["force"]
        else:
            raise ValueError("No deflection data found!")

        if "data" in a:
            b =  a["data"]
        elif "deflection" in a:
            b = a["deflection"]
        elif "raw_data" in a:
            b =  a["raw_data"]
        else:
            raise ValueError( "Bogous force data!" )

        if self._indx_.size > 0:
            return b[self._indx_]
        else:
            return b
    #end force

    @property
    def force_unit(self):
        """ Get the default force unit stored in the headers"""
        a = self.full_segment
        if 'force' in a:
            a = a['force']
        if "unit" in a:
            return a["unit"]
        elif "deflection_unit" in a:
            return a["deflection_unit"]
        elif "raw_unit" in a:
            return a["raw_unit"]
        else:
            return ""
    #end force unit

    @property
    def force_constant(self):
        """ get or set the force constant for calculating forces from
            deflection
        """
        a = self.dataset

        if 'force constant' in a:
            return a['force constant']
        else:
            return -1.0
    #end get force_constant

    @force_constant.setter
    def force_constant(self, k):
        a = self.full_segment

        if 'force' in a :
            a = a['force']

            if "deflection" in a:
                a['data'] = a['deflection'] * k
            else:
                raise ValueError("No deflection is defined!")

            #force constant is not inside the segment:
            self.dataset['force constant'] = k

            if 'deflection_unit' in a:
                if a['deflection_unit'] == 'm':
                    a['unit'] = 'N'
                elif a['deflection_unit'] == 'nm':
                    a['unit'] = 'nN'
                else:
                    a['unit'] = "%sN" %(a['deflection_unit'].split('N',1)[0])

        else:
            raise ValueError("There is no force data available!")
    #end force_constant

    @property
    def Zmax(self):
        """ maximal limiting Z value to be used
            Setting this to > 0.0 will cause data being filtered
            to (Z - Z.min()) <= Zmax
            setting <= 0 will turn the filter off
        """
        return self._Zmax_
    #end Zmax general

    @Zmax.setter
    def Zmax(self, zm):
        #reset before applying the filter:
        self._indx_ = nu.empty(0)
        self._Zmax_ = zm
        #keep Zmax, because it may be meaningful for other segments!

        # short process: if negative or zero
        # it is now removed, so return:
        if zm <= 0.0:
            self._Zmax_ = -1.0
            return

        Z = self.Z
        Z = Z - Z.min()

        if zm > Z.max():
            #nothing to do
            return
        #now, we do sg.:
        self._indx_ = Z <= zm
        return
    #end setter Zmax

    @property
    def Z(self):
        """ Piezo distance measured or direct if measured is
            not available
            The local value of Z0 is automatically subtracted,
            so you can use it to define the contact point.
        """
        a = self.full_segment
        if "height" in a:
            a = a["height"]
        else:
            raise ValueError("No Z-piezo information found")

        if "data" in a:
            Z = a["data"]
        else:
            raise ValueError("Invalid piezo data set")

        if self.Z0 is not None:
            Z = Z - self.Z0

        #return the result (filtered if index is set)
        if self._indx_.size >0:
            return Z[self._indx_]
        else:
            return Z
    #end of Z (piezo height)

    @property
    def Z_type(self):
        """ is this a measured height or a piezo height? """
        a = self.full_segment
        if 'height' in a:
            a = a['height']
        else:
            return ""
        if 'type' in a:
            return a['type']
        else:
            return ""
    #end Z_type

    @property
    def Z_unit(self):
        """ unit of piezo position from the segment data """
        a = self.full_segment
        if "height" in a:
            a = a['height']
        else:
            return ""
        if "unit" in a:
            return a["unit"]
        else:
            return ""
    #end Z_unit

    @property
    def t(self):
        """ time axis for the data
        """
        a = self.full_segment

        #this variable is not clear yet, I need test data on this!
        t = nu.empty(0)
        if 'time' in a and "data" in a['time']:
            t = a['time']['data']
            if self._indx_.size > 0:
                t = t[self._indx_]

        return t
    #end t

    @property
    def sensor_response(self):
        """ find and return the sensor response
            In setting this property the deflection and force are
            recalculated automatically.
        """
        a = self.dataset
        if "sensor response" in a:
            return a['sensor response']
        else:
            return -1.0

    @sensor_response.setter
    def sensor_response(self, sn):
        """ Set the sensor response and recalculate the
            deflection and the force data
            Unit: the distance unit after raw_data*sensor_response
        """
        a = self.full_segment
        if 'force' in a:
            a = a['force']
        else:
            raise ValueError("Can not set sensor response with no defleciton!")

        #if there was no sensor response defined, we do that now:
        if sn > 0.0:
            self.dataset['sensor response'] = sn
        else :
            raise ValueError('Invalid sensor response %.3f' %sn)

        if 'raw_data' in a:
            a['deflection'] = a['raw_data'] * sn
            if "deflection_unit" not in a:
                if nu.log10(sn) < -5:
                    a["deflection_unit"]= "m"
                elif nu.log10(sn) < -2:
                    a['deflection_unit'] = 'microm'
                else:
                    a['deflection_unit'] = 'nm'

            k = self.force_constant
            if k > 0.0:
                a['data'] = k * a['deflection']
        else:
            raise ValueError("Invalid deflection data, no raw signal found!")
    #end sensor_response

    @property
    def deflection(self):
        """ Deflection data.
            You should not set deflection directly, it is derived
            from the raw data by setting the sensor_response!
        """
        a = self.full_segment
        if "force" in a and 'deflection' in a['force']:
            defl  = a["force"]['deflection']
        else:
            return nu.empty(0) 

        if self._indx_.size > 0:
            return defl[self._indx_]
        else:
            return defl
    #end deflection

    @property
    def deflection_unit(self):
        """What unit the deflection has """
        a = self.full_segment
        if "force" in a:
            a = a['force']
        else:
            return ""

        if "deflection_unit" in a:
            return a['deflection_unit']
        else:
            return ""
    #end deflection_unit

    @property
    def deflection_raw(self):
        """ The measured deflection data, in V.
            This is the basis of calculating the deflection in
            length (distance) units as well as the force.
            Thus, if you modify this, all changes.
            Here you can apply base line correction and filtering.
        """
        a = self.full_segment
        if "force" in a:
            a = a["force"]
        else:
            print("vertical deflection was not measured!")
            return nu.empty(0)

        #now we are within data['vDeflection'][...]
        if 'raw_data' in a:
            if self._indx_.size > 0:
                return a['raw_data'][self._indx_]
            else:
                return a['raw_data']

    @deflection_raw.setter
    def deflection_raw(self, data):
        """ you can change the raw data if needed...
            (for the sake of completeness)
            Specify a full set, not only the filtered one!
        """
        a = self.full_segment

        #turn off filtering for a moment:
        ZmaxB = self.Zmax
        self.Zmax = -1.0

        #set, if it fits to the distance array:
        if 'force' in a \
                and data.size > 0 \
                and len(data) == len(self.Z):
            a = a['force']
            a['raw_data'] = data

            #now we have to propagate the changes:
            sn = self.sensor_response
            if sn > 0:
                #the setter generates the deflection data,
                #direct setting is not allowed
                #a.sensor_response = sn
                #alternatively we could write directly (internal use):
                a['deflection'] = sn*data

                k = self.force_constant

                if k > 0:
                    a['data'] = k* a['deflection']
        #end if vDeflection
        self.Zmax = ZmaxB
    #end deflection_raw


    @property
    def distance(self):
        """ The tip distance from the surface.
            This is only a relative value if the surface location
            and the base line are not corrected!
            (It is a shortcut for Z - deflection)
        """
        return self.Z + self.deflection
    #end of distance

    def invert(self):
        """ Invert the sequence of all data arrays """
        a = self.full_segment

        #inverting may mess with the Zmax settings
        #thus, reset it:
        zmax = self._Zmax_
        self.Zmax = -1 #this takes care of turning it off

        z = self.Z
        # check order of the array, if it was not empty
        if z.size > 0:
            if (z[-1] - z[0]) < 0.0:
                print("needs to be inverted")
                #run through the segment keys, and find data within
                #targets are: ['force']['deflection'], ['height']['data'], etc.
                for i in a.keys():
                    b = a[i]
                    if type(b).__name__ == 'dict':
                        for j in ['raw_data', 'data', 'deflection']:
                            if j in b :
                                b[j] = b[j][::-1]
                    #end if dict
                #end for in keys
            #end if invert
        #end if no z

        #reapply:
        self.Zmax = zmax
    #end invert

    def baseline(self,
                 end=0.05,
                 Zrange= -1.0,
                 order= 0,
                 endN = -1,
                 wave = False,
                 Nsmooth= 30,
                 verbose=False
                 ):
        """ Calculate a baseline using the Baseline() generic
            function.
            Parameters:
            end         what relative length of the curve should be
                        used at the end
            Zrange      the same in length, use the same units as Z
                        this redefines endN if > 0.0

            order       the order of polynomial to use. 0 stands for
                        average, 1 for linear etc.
            endN        an alternative to end, specify the number of
                        points instead of portion
                        This takes over if > 2
            wave        Boolean, if True, use a wave background not a
                        polynomial one. If set, order is not used.
            Nsmooth     number of points in the smoothing filter of wave bg.
            verbose:    give feedback

            We apply Baseline to the Zmax limited extension, then
            recalculate the full background corrected deflection_raw
            array!

            For further details see the Baseline() function.

            return:
            returns the base dict and corrects deflection_raw,
            deflection and force accordingly.
        """
        #here the Zmax takes effect in self.Z...
        Z = self.Z

        if Zrange > 0.0 and endN < 2:
            z = self.Z #thus Zmax and Z0 are applied

            if z.min() > z.max() - Zrange:
                raise ValueError("Too broad range, curve is only %.2f long"\
                            %(z.max() - z.min()))
            #end if too large
            #convert zmax to an index:
            endN = int((z > (z.max() - Zrange)).sum())
            print("New endN is %d" %endN)
        #end if end is invalid

        if wave:
            N = len(self.Z)
            i1 = N - endN if endN > 2 else int(N*(1.0 - end))
            base = fit_wave_background(self.Z[i1:], self.deflection_raw[i1:], Nsmooth)

        else:
            #call the generic baseline correction:
            base = Baseline(self.Z, self.deflection_raw, order= order, \
                    end=end, verbose=verbose, endN = endN)

        #turn off the filter:
        ZM = self.Zmax
        self.Zmax = -1
        if wave:
            bg= bg_wave(self.Z,
                        base['amplitude'],
                        base['amplitude_slope'],
                        base['omega'],
                        base['beta'],
                        base['delta'],
                        base['offset'],
                        base['slope']
                        )
        else:
            bg = nu.polyval(base['fit'], self.Z)

        self.deflection_raw = self.deflection_raw - bg
        sens = self.sensor_response
        if sens > 0:
            self.sensor_response = sens #recalculate deflection

            k = self.force_constant
            if k > 0:
                self.force_constant = k #this will recalculate force

        #restore the filter:
        self.Zmax = ZM
        return base
    #end of baseline

    def find_contact_point(self):
        """ identify the contact point for this segment using the
            Contact_point generic function.
            Sets Z0, which is subtracted from Z automatically
        """
        self.Z0 = Contact_point(self.Z, self.deflection_raw)
    #### CONTINUE

#end class AFMforce

###############################################################################
########## General functions for conversion and analysis
###############################################################################


def Baseline(z,I, end=0.05, order=0, verbose=False, endN=-1):
    """ Do a baseline correction on the intensity data.
        Ideally every measured force curve ends in a constant, zero
        force tail at high distances. If the measurement did not reach
        this region, this function will generate garbage.
        In some occasions reflections or other sources may result in a
        linear background, which would make the identification of adhesion
        forces difficult.

        The function will first smoothen the force curve, then calculate
        a difference data set (f[1:] - f[:-1]). Then use the end part
        (defined by the end variable (0-1)), and estimate the standard
        deviation of this part. Then we seek the first point from the
        back, where the difference curve reaches more than 3 times this
        standard deviation. We can apply a polynomial fit from this point
        to the end. This is the background being subtracted.

        Parameters:
        z           position data
        I           vertical deflection data
        end         how much of the dataset to be considered as the end
                    part for the standard deviation. It should not be
                    more than about 10%!
        order       the order of the polynomial fitted. If 0, then
                    the average is taken. It should not be more than 1 in
                    general, but 0 is also sufficient in common cases.
        endN        specify end in number of points
                    if > 2, it takes preference instead of end.

        return:
        a dict containing:
        'data'      the corrected I data
        'z'         a copy of the z data
        'bkg'       the evaluated background values
        'fit'       the fit array from numpy.polyfit or the average
    """

    if z.size != I.size:
        print("non-matching z,I data!")
        return {}

    if z.size < 100:
        print("The data set may be too short!")

    #we allow the user provide number of points instead of portion:
    Nend = int(end*z.size) if endN < 2 else int(endN)
    #alternatively use the Z values for defining Nend:
    #Nend = (z < end*z.max()).sum() if endN < 2 else endN

    if Nend < 5:
        print("Warning, low number of end points %d" %Nend)

    if order > 1:
        print("Warning: more than linear polynomials may falsify your data!")
    #end ifs

    #check for order (we assume no jump back btw. start and end:
    if z[-1] < z[0]:
        z = z[::-1]
        I = I[::-1]
        if verbose:
            print("swapped order!")

    # this part was ineffective.... (not used)
    #
    #if Nend < I.size:
    #    I0 = I[-Nend:].mean()
    #    ff = I - I0
    #else:
    # and this should not happen
    if I.size >= Nend:
        print("Invalid number of end part! %d" %Nend)
        return {}
    #end if

    #do the fit to the end:
    ff = nu.polyfit(z[-Nend:], I[-Nend:], order)
    #print ff.shape

    # expand the fit as close to zero as it gets
    # to do this, we investigate how the noise
    # changes

    #evaluate the fitted curve and see the error:
    df = (nu.polyval(ff, z) - I)**2
    #the average squared error at the end:
    sigma = df[-Nend:].mean()

    # if the end is noisy, we seek the first point
    # before that part
    indx = (df[:-Nend] > 3*sigma).nonzero()[0]
    if len(indx) < 1:
        print("no hits found, all df is constant!")
        i0 = 0
    else:
        i0 = indx[-1]
    #end if
    ########## This may be too much here....
    #if the result is too close to the start, it may be completely off
    #fall back to the standard way:
    if i0 < 0.2*len(I):
        print("suspicous low i0: %d" %i0)
        i0 = I.size - Nend
    #end i0

    if order > 0:
        fit = nu.polyfit(z[i0:],I[i0:],int(order))
        fitted = nu.polyval(fit, z)
    else:
        #fit shold be an array, even if it is 1 long
        fit = nu.asarray([I[i0:].mean()])
        fitted = nu.zeros(I.size)+ fit[0]
    #end if order

    if verbose:
        fig = pl.figure(1)
        # go for the wrappers and let matplotlib
        # do its job
        pl.clf()
        pl.plot(z,I, 'b+', alpha=0.8)
        pl.plot(z,fitted, 'g-+', alpha=0.3)
        pl.xlabel('Z position')
        pl.ylabel('vertical deflection')
        pl.draw()
        print("Fit parameters:")
        print("order: %d" %order)
        print("fit:", fit)
        print("var: %.4f, Nend: %d" %(sigma, Nend))
        print("first point: %d, %.3f" %(i0, z[i0]))
        print("number of points: %d" %I.size)

    return({'data':I-fitted, 'z':z, 'bkg':fitted, 'fit':fit,
            'i0':i0, 'sigma': sigma})
#end Baseline


def Contact_point(z, f,
            radius= 10,
            par= 2,
            method= 'SavGol',
            use_deriv= True,
            noiselevel= 0,
            verbose=False):
    """ Identify the contact point from a piezo distance vs.
        deflection data. Best is to call it with the raw sensor
        data, but it is up to the user.
        However, having the baseline correction may improve
        the result.

        If derivative is used, the algorithm takes a smoothened deflection,
        then define the contact point as the smallest piezo position
        where the deflection gradient is non-negative
        Close to zero has to be allowed for curves without adhesion
        (negative) part.

        Otherwise take the lowest Z value where the curve is
        negative or below noiselevel
        (use baseline correction before calling this!)
        **Use a linear interpolation to improve the Z estimate.**

        This should work in most cases, but be careful with data
        not hitting the wall...

        Parameters
        z           distance or Z-piezo data
        f           force or deflection values
        radius,
        par,
        method      are passed to Smooth
                    if radius is <= 0, do not use smoothing
        use_deriv   use a derivative rule above, or the zero crossing (false)

        noiselevel  if zero crossing is used, this is the level below
                    zero, where we consider the curve negative.
                    It is possible that noiselevel > 0, then it is the leftmost
                    point where the curve has deviated this much from zero.
                    Sometimes this is used to estimate the point of contact
                    for indentation.


        verbose     plot the result indicating the smooth curve and the
                    estimated contact point

        return:
        the z0 value.
    """
    N = len(z)

    if N != len(f) or N < 2:
        raise ValueError("Invalid length of z and/or f data!")
    # end if z length
    if z[0] > z[-1]:
        print('inverting sequences!')
        z = z[::-1]
        f = f[::-1]

    if radius >0:
        fs = Smooth(f, radius= radius, par= par, type = method)
    else:
        fs = f

    if use_deriv == True:
        # difference of the smoothened curve should be negative
        # as we are heading from the wall to the end
        # loading the force ensures the order: z[-1] > z[0]

        # df = fs[1:] - fs[:-1]
        if radius > 0:
            df = Smooth(f, radius= radius, par= par, type = method, deriv= 1)
        else:
            df = f[1:] - f[:-1]

        indx = (df >= 0).nonzero()[0]
        # print(indx)
        if len(indx) < 1:
            print('Warning: no contact point is found')
            return 0.0
    else:
        indx = (fs < noiselevel).nonzero()[0]
        if len(indx) < 1:
            print('Contact point not found!')
            return 0.0

    # both methods generated an index, from which we get:
    z0 = z[indx].min()

    if verbose:
        pl.clf();
        pl.plot(z - z0,f,'b+')
        pl.plot(z - z0, fs, 'r-')
        pl.plot( (0,0), (f.min(), f.max()), 'g-')
        pl.draw()
        print("contact point found as: %f" %z0)

    return z0
#end Contact_point

def Find_puncture(z, f,
                 radius= 20.0, width= 5.0, method='Gauss',
                 verbose = False):
    """ use the derivative (index) of the force curve to estimate
        jump regions. It works the same way as the derivative in contact
        point estimation above, but looks for the various segments
        with the same property: where the derivative becomes non-negative.

        Since the smoothing will shift the position, after finding the
        candidate, take the local maximum of the original data, +/- 2 width
        points around the detected candidate point.
        If multiple maxima occur, take the first one.

        To estimate the jump, search for a local minimum to the left.
        Use 2x width first, but shift on if the minimum is at the left
        edge of the window. This method assumes the force curve is increasing
        to the left as it is hitting the wall, except for local minima at
        rupture points or jump into contact.

        Parameters:
        z,f         a force data set, possibly sorted
        radius      window radius of the smoothing curve
        width       width of smoothing kernel; both passed to Smooth
        verbose     provide some plots

        return:
        a dict containing:
        'indx':     index of the positions in the provided array
        'z'         corresponding z values
        'f'         corresponding force value
    """
    indx = (Smooth(f, radius, width, type= method, deriv=1) >= 0).astype(int)
    if indx.sum() < 1:
        print('Only negative slopes were found')
        return nu.ndarray()

    fs = Smooth(f, radius, width, type= method)

    indx2 = indx[1:] - indx[:-1]

    res = (indx2 < 0).nonzero()[0]
    N = len(f)
    width = int(width)

    newres = nu.zeros(res.size, dtype= int)
    minres = nu.zeros(res.size, dtype= int)

    for i in range(res.size):
        ii = res[i]
        i0 = max(0, ii - 2*width)
        i1 = min(N, ii + 2*width+1)
        newres[i] = (f[i0:i1].max() == f[i0:i1]).nonzero()[0][0] + i0 if i0 != i1 else i0

        i0 = max(0, ii - 2*width)
        ii = (fs[i0:ii]  == fs[i0:ii].min()).nonzero()[0][0] + i0 if i0 != ii else i0

        # if the minimum is further off to the left, follow it:
        # the curve should be increasing in this direction in general
        # search in the smooth curve, to avoid noise resulted minima
        while (ii == i0 and ii > 0):
            i0 = max(0, ii - width)
            ii = (fs[i0:ii]  == fs[i0:ii].min()).nonzero()[0][0] + i0

        # and refine locally:
        i1 = min(N, ii + 2*width+1)
        minres[i] = (f[i0:i1] == f[i0:i1].min()).nonzero()[0][0] + i0

    if verbose:
        pl.clf()
        pl.plot(z, f, '-+')
        pl.plot(z, fs, 'r-', alpha= 0.5)
        pl.plot(z[newres], f[newres], 'ro')
        pl.draw()

    return {'indx': newres, 'z': z[newres], 'f': f[newres],
            'dz': z[newres]-z[minres], 'df': f[newres]-f[minres]}
# end of Find_puncture

def Z_to_D(z, f, factor=1.0, verbose=False):
    """ Do the conversion between raw data and force-distance curve.
        This is a basic step for force measurements.
        Separating this step from identifying the parameters allows
        the user to decide what parameters to use. For soft surfaces,
        or curves with too short contact or zero force parts, the
        user may supply data from other measurements.

        ** this is more a legacy function to be used when the force class
           is not suitable for some reason **

        Parameters:
        z           height data (piezo or strain gauge...)
        f           vertical deflection (sensor volts or force values)

        factor      multiplier constant: (watch the dimensions!)
                    thus f = z*factor for the hard wall
                    it should be >0, if not, we make it so
                    (it is the same here if it is N/m or V/m... etc.)

        Return
        the corrected distance values: d = z + deflection

    """
    if len(z) != len(f):
        raise ValueError("The length of input arrays does not match")
    #end if

    if factor < 0.0:
        factor = -1.0*factor
    #end if

    #the actual distance of the tip (probe) from the surface is then:
    defl = f/factor + z if factor != 0 else z

    if verbose:
        fig = pl.figure(1)
        fig.clf()
        # from 2020 go for the wrappers
        #plt = fig.add_subplot(111)
        #plt.plot(z, f, 'r+', alpha=0.8)
        #plt.plot(defl,f,'b+', alpha=0.3)
        #plt.plot(z, nu.zeros(defl.size), 'g-', alpha=0.3)
        #plt.set_xlabel('Z position')
        #plt.set_ylabel('vertical deflection')

        pl.plot(z, f, 'r+', alpha=0.8)
        pl.plot(defl,f,'b+', alpha=0.3)
        pl.plot(z, nu.zeros(defl.size), 'g-', alpha=0.3)
        pl.xlabel('Z position')
        pl.ylabel('vertical deflection')

    return defl
#end Z_to_D

def find_slope(z,f, factor=3.0, Nfit= 5, verbose=False):
    """ Find the initial slope of the force curve, based on a simple
        linear fit. The data should be base line corrected first!

        Fit the first segment in Nfit subsegments, and find the
        most negative slope of them. Then extend the fit to the
        maximal range and return the result.
        The error of the fits is analyzed, if they show a large
        variation, then ignore the most noisy 1/4 fits, and use
        only the rest.
        If the data is consistenly noisy, this does not change a thing.

        Parameters:
        z, f:        the dataset: z-piezo and force (raw) data
        factor:     a multiplier of the error to estimate
                    deviation from the base line. (3 or 5 are suitable)
        Nfit:       number of subsegments
        verbose     provide some plot

        return
        the absolute value of the slope (which is <0 originally)
    """
    #First we start finding a threshold: where do we reach
    #the zero force line -> the slope may extend up to that point
    Nend = int( max(0.05*len(f),10))

    if len(f) < 5*Nfit:
        print("too short data set %d" %len(f))
        print("or set Nfit to: %s" %(Nfit/5))
        return 0.0
    #end if

    #if one wants stupid results, let them have fun...
    if factor < 0 :
        factor = -1.0 * factor

    sigma = f[-Nend:].std()
    if nu.abs(f[-Nend:].mean()) > 0.1:
        print("warning, data is probably not base line corrected!")
    #end if

    #Now we can estimate which part to look at: 0 - i:
    indx = (f < factor*sigma).nonzero()[0]
    if len(indx) < 1:
        print("no crossing found, terminating")
        return 0.0

    iend = indx.min()
    if iend < 3*Nfit:
        print("Insufficient data length: %d vs. %d" %(iend, 3*Nfit))
        return 0.0
    #end if

    #now fit segments:
    fits = []
    errs = []
    di = int(iend/Nfit)

    for i in range(Nfit):
        i0 = i*di
        i1 = i0 + di
        fit = nu.polyfit(z[i0:i1],f[i0:i1],1)
        if fit[0] > 0:
            # positive slope means we have something weird, not a
            #sensor response part
            #print("invalid fit %d, skipping" %i)
            continue
        else:
            fits.append(fit)
        err = ((f[i0:i1] - nu.polyval(fit, z[i0:i1]))**2).mean()
        errs.append(err)
        #print i,":",err
    #end if

    if len(errs) < 1:
        print("Invalid data set")
        return 0.0

    #estimate a maximal error level:
    # errm * factor
    errm = min(errs)

    #is the noise inhomogeneous?
    # the noise range should be less than factor*errm
    if (max(errs) - errm) > factor*errm :
        # we sort the errors, but also the corresponding fits:
        err_sorted = sorted(enumerate(errs), key= lambda x: x[1])
        errs = [i[1] for i in err_sorted]
        fits = [fits[i[0]] for i in err_sorted]

        # take the best 75% of the list:
        end_indx = int(0.75*len(errs))
        limit = errs[end_indx]
        # collect those who are below this
        if end_indx > 0:
            errs = errs[0:end_indx]
            fits = fits[0:end_indx]
    #end if

    # get the most negative slope (steepest, but negative):
    fit = min(fits, key=(lambda x: x[0]))
    # the elements are numarray objects, index will not work
    # find the one which was the minimum:
    i = [i for i, f in enumerate(fits) if f[0] == fit[0] and f[1] == fit[1]]
    i = i[0]

    errm = errs[i]

    # reevaluate the steepest fit to find where it matches the best:
    err = (f[:iend] - nu.polyval(fit, z[:iend]))**2
    indx = (err <= errm).nonzero()[0]
    #if it is valid, refit:
    if len(z[indx]) > 2:
        fit = nu.polyfit(z[:iend][indx],f[:iend][indx],1)
    #else we keep the minimum, without refit

    else:
        print("refit failed with index:", indx)

    if verbose:
        fig = pl.figure(1)
        fig.clf()
        #plt = fig.add_subplot(111)
        #plt.plot(z,f,'r+')
        #plt.plot(z[:iend],nu.polyval(fit,z[:iend]),'g-')
        pl.plot(z,f,'r+')
        pl.plot(z[:iend],nu.polyval(fit, z[:iend]),'g-')
        pl.draw()

        print("slope finding parameters")
        print(f'Start with index: {i} fit: {fits[i]}, error: {errs[i]}')
        if indx.size > 2:
            print(f'fit index range: {min(indx)}, {max(indx)}')
        else:
            print('reevaluation came up empty')
        print("fit:", fit)
        print("in:", fits)
        print("belonging error:", errm)

    return nu.abs(fit[0])
#end find_slope


def FindMinima(f, dN= 50, Z= nu.empty(0), dZ = -1.0, minf = None, verbose=False):
    """ find the local minima in a force curve
        Use a simple local minimum search, then check that
        this value is minimal in the next dN points as well.
        Thus, if the curve ends with a minimal plateau, it may
        miss the last point.

        Parameters:
        f       force data (ndarray)
        dN      the number of data points to consider for a local
                segment
                If dZ is specified, it is taken priority. However, if
                dZ results in < 3 points, dN overrides with a warning.
                This is to prevent premature stop of the search due to
                too few data points.

        Z       the distance data
        dZ      window size along the distance data
                if > 0 and Z is defined, it takes preference
        minf    if defined, only only f values below this are accepted

        verbose     show some feedback

        Return value:
        a list of indices indicating where the local
        minima are found
    """
    N = len(f)
    dN = int(dN)

    if N < 1 :
        raise ValueError("Empty parameters!")
    if len(f) < dN:
        return (f == f.min()).nonzero()[0]
    #end ifs

    if len(Z) > 0 and len(Z) == N:
        UseZ = True
    else:
        UseZ = False
    #end if len(Z)

    sum = nu.sum

    i=0
    currpos = 0
    minlist = []
    # hit marks if we have found a local minimum:
    # Going from left to right,
    # if the first minimum is on the left edge, it is a local one
    hit = True
    while (i < N):

        #set up the range of interest:
        if UseZ:
            Z1 = Z[i] + dZ
            i1 = sum(Z < Z1)
            #problem: if dZ is too small, this stops the loop:
            if i1 < (i+3):
                print("Warning, too short dZ, overriding at %d!" %i1)
                i1 = min(i+dN, N)

        else:
            i1 = min(i+dN, N)
        #end if UseZ

        #do we still have data to the right:
        if (i1-i) < 3:
            break
        #end if end of data

        #we find the minimal points in the segment:
        indx = (f[i:i1] == f[i:i1].min()).nonzero()[0]
        #and take the rightmost one:
        imin = indx[-1]
        #now, we have a local minimum
        #is it a new one, or a previously stored one?
        #starting from the beginning, we have one at the first
        #point or somewhere between
        #if we put that to i, the next segment may have a lower
        #point or not...
        #One thing to check: is the curve flat?
        #   THen what is the minimum?
        if imin > 0:
            #we have a minimum further up, it must be real:
            hit = True
            i += imin
            #temporary storage of the location
            currpos = i
        #if not, then it is the leftmost point.
        elif hit:
            #this is a real hit, so store it from the temp storage:
            if minf is None or (f[i:i1].mean()-f[currpos]) > minf:
                minlist.append(currpos)
                #if it is a hit, then reset the search:
                hit = False

            #jump to the end, and go on searching:
            i = i1
        else:
            #it is no hit, but i=0
            #jump through, this is an upwards slope:
            i = i1
        #end if
    #end while i<N
    #check for the last point:
    if imin+i == N-1:
        minlist.append(N-1)

    if verbose:
        print("found %d local minima" %len(minlist))
        fig = pl.figure(1)
        fig.clf()
        #plt = fig.add_subplot(111)

        if UseZ:
            pl.plot(Z,f,'+')
        else:
            pl.plot(f, '+')
        fmi = f.min()
        fma = f.max()
        for i in minlist:
            pl.plot([i,i],[fmi,fma], 'g-')

        pl.draw()
        print("Last range: %d:%d" %(i, i1))
        print("Last minimum: %d" %imin)
    #end verbose
    return minlist
#end of FindMinima


def Indentation(Z, Defl, z0=0, deflection0=0, defN= 50,\
        auto=False, verbose=False):
    """ Calculate the indentation from the deflection and the z position
        data of a segment.
        Let deflection0 be the background constant deflection,
        Z0 the Z (piezo) position of the contact point.

        The indentation (D) is: Z0 - Z - (deflection - deflection0) for
            Z  < Z0, and 0 for all Z >= Z0.

        This can be converted to: D0 - Z - deflection, where
            D0 = Z0 + deflection0

        Thus, if we do not know Z0 and/or deflection0, we have only
        a D0 offset. While normally the D < 0 has no meaning, and invalid,
        we keep those values, so the user may correct for various D0 values
        e.g. in fitting elastic models to the force curve.

        Automatic: the standard deviation of the deflection is calculated
        (noise level) of the last 10% of the curve but at least 50 poitns.
        Then the lowest Z value where the deflection > negative standard
        deviation is the approximated Z0 point.
        Rationale: the curve starts deviating, then the noise may still cross
        0 up and down. This should be more sensitive than picking the
        Z-position of the lowers zero crossing.

        Parameters:
        Z           the Z piezo data
        Defl        the deflection data (distance)
        z0          the zero indentation position on Z
        deflection0 see above

        defN        how many points to use in the estimation of the
                    zero defleciton?
        auto        True or False
                    if True, try estimating a Z0 value and a deflection0

        return:
        a dict containing:
        "D"         the deflection values
        "D0"        the offset
        "Z0"        the used Z0
        "def0"      the used deflection0

    """
    if (Z[-1] - Z[0]) <  0:
        Z = Z[::-1]
        Defl = Defl[::-1]

        if verbose:
            print("swapped order!")
    #end if swap order

    if auto:
        Nauto = int( max(len(Defl)/10, defN) )
        defl0 = Defl[-Nauto:].mean()
        th = Defl[-Nauto:].std()
        if verbose:
            print("zero deflection", defl0, "threshold:", -th)

        # the first point where the noisy data
        # goes below 0
        z0 = Z[ (Defl-defl0) < -th ].min()
    else:
        defl0 = deflection0
    #end if

    indentation= { "D": (z0 + defl0 - Z - Defl), "z0":z0, \
            "D0": z0+defl0, "def0": defl0, "z": Z - z0}
    if verbose:
        print("Estimated Z0:", z0)
        print("Deflection 0:", defl0)

    return indentation
#end of Indentation

