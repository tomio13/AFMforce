Disclaimer:
    This program has no connection to the JPK Instruments AG or Bruker.


Test files are from:
/home/tamas/Projects/AFM/Integrin-ligand-interaction/2012-07-19/1st/
time stamp: 22:18:56.
License: Creative Commons (CC BY 4) https://creativecommons.org/licenses/by/4.0/

2012-09:
    I need algorithms to analyze our molecular recognition force data. The
    task is multifold:
     - import JPK force data possibly from the raw files
     - extract calibration information
     - convert the data to force-distance curves
     - identify adhesion segments
     - extract maximum adhesion forces
     - extract location of the experiments
     - extract bond strenght from jumps
     - fit worm-like-chain model if possible
     - do batch processing and informative output

    The calibration done in air is valid to the spring constant, but the block
    of liquid may change the sensore response slightly. If the data is good
    enough, it is better extracting the sensor response from the actual raw
    data.

jpk-force:
    The raw file is a zip file containing subfolders, text headers and binary
    dumps of the raw data. Analyzing the text headers show that all contain
    various conversion information. The file structure seems quite
    straightforward. It took a day or two, but the script can nicely read the
    data now and dump it into a dict structure.

2013-03-22:
    Well, clean up the code a bit, and figure out how to get short integer
    types to work.

    The importer is in alpha phase now. It reads much more of the header and
    figures out some of the correlated parameters into dicts. The ultimative
    goal is to properly parse the header and use all comments on what is
    related to what.
    Better explaining by example: The experiment type is specified in the
    first part of header.properties under the type variable. Right now this is
    not used, though one can pick all further experiment parameters using the
    value of type. Thus, the corresponding description and header are found
    below that line.

    Now in the segment-header.properties it is used however, to group the
    calibration parameters as dicts to the specific segments/channels.
    One has to pick segment-0...segment-3 from the data dict, then for example
    height, and there one can find what conversions are available with what
    parameters, together with the raw data (key: data).

    Thus, a channel measured in volts may be converted to nominal distance,
    absolute distance etc. if these are specified in the header.

2013-04-17:
    Modified the slope finding routine to a simpler approach: chop up the
    initial segment to N parts, each > 5 pixels long. (N is user defined.)
    Then fit them, and use the steepest negative slope as the hit.
    Refine the fit then so, that the most nearby points are included, by
    analyzing the error of the data in the whole first segment.

    If the error of the fitted parts varies a lot (factor*min(errors)), then
    drop the worst 1/4th to improve the result. This may not be optimal, just
    a start to decrease confusion caused by high contact noise.

2013-11-20:
    The next version, 2.0 contains a slightly different header structure. I
    have rewritten the whole import system to play nicer both with the old and
    new file format.
    The key point: header_to_dict() now makes an embedded dict structure,
    converting the a.b.c headers to a{b:{c:value}} dicts. This results in a
    well structured header.
    The CalibrateForce() is also adjusted to handle this. The only critical
    part is to chase down the encoder parameters. Maybe this should be moved
    to the reading routine, I am not sure at the moment.
    Anyway, the code has tested out fine for the old and new force data.

2014-04-09:
    Start working on indentation data. First, calculate indentation from the
    deflection and Z-position. A simple loop can estimate a zero indentation
    distance, but this is somewhat crude. It may be necessary to let models
    decide this later.
    Anyway, the Indentation() function allows doing simply the subtraction, or
    can handle some estimations itself.

    The first cone model fitting is ready to trial.

2014-10-28:
    a bit of bug fixing: now the fit_Hertz should work fine for spheres too

2014-12-10:
    Fixed the fitting somewhat. There are two things fixed:
    - when calling the fit, it is better to have the zero Z a little deeper,
      so the fit starts in a broader range
    - the second the signs of the guessing is fixed, now if we take an
      absolute value of the first parameter, then we do that in every use of
      it. THis has fixed some occasional missing of the curves. The algorithm
      now does not throw a singular value fit in every 10 cases.

2015-03-09:
    Stavisky-Golay filter to improve data smoothing. This one is a strong
    tool, which can keep the position of the peaks. Important for adhesion
    data smoothening. It is now part of the Smooth command, and the default
    filter form.
    Parameters: a half length of the smoothening window (radius), and the
    order of interpolation (1 o3 2 are suitable). It is more versatile than as
    it is used here though with possible derivatives to be calculated.

2015-03-26:
    Added a WLC fitting algorithm. This is the simplest WLC model to estimate
    a contour length and a correlation length.

2015-05-20:
    Start up with a JPKforce class. This is to simplify all the data treatment
    doing automatically the conversions between raw data, force, indentation,
    Z-position and distance, allowing for base line correction and finding the
    approximate contact point.
    It uses plenty of properties to do this all quite automatically, removing
    many of the complexity from the processing script.

    Added a Contact_point() general function to estimate the contact point
    based on the difference in the smoothed deflection data.
    It will not work on curves with punctures or weird increasing contact
    parts (increasing with distance, not decreasing as usual).
    Setting Z0 manually is also possible.
    In the class, one can not set force and deflection manually, only the
    deflection_raw data, which is the sensor values in volts. Changing the
    sensor_response and / or spring_constant will recalculate the deflection
    and the force values automatically.
    Selecting a new segment takes care of setting the data in increasing
    Z-piezo series.

    The original data is still available in dicts inside the class, such as:
    JPKforce.raw_dataset and JPKforce.calibrated_data
    JPKforce.full_segment is the dict of the actual segment.

2015-05-27:
    added a time variable, to extract time dependence. This way we should be
    able to estimate a local pulling rate from the WLC fits as well as from
    the data.

    Dropped the F0 parameter from the WLC fits. It should not be there,
    and fitting it makes more trouble than it worths. In some cases the fit
    will look quite bad because of this, but in general the result is more
    reliable (physically).

2015-06-29:
    Minor corrections to the class to clean the code up a bit.
    The FindMinima() is extended. Now it allows for providing a Z array and a
    dZ window size (half of it), which it uses then to identify the
    corresponding dN values. This makes it a bit more dynamic, because at
    pulling forces there are more points in the same dZ piezo range than in
    relaxed parts.
    
    The Baseline is also accepting direct dN value instead of the percentage
    only. Since in the correction the piezo movement is expected to be
    constant velocity, so the dN should directly correspond to a dZ here.
    (I may implement a dZ value as well.)

2016-01-24:
    Minor bug: numpy 1.10 does not recognize i32 anymore, it is i4 (byte size,
    not bit size). 

2016-09-23:
    A bug in the force constant setting part. It did not write the right field
    in setting a new constant. And it has to set the unit, because it may be 
    still in deflection units.

2017-03-20:
    There was a bug in the SavGol filter I did not notice before, but was not
    affecting the performance of the code. The result has a multiplier:
    (-1)**difforder, thus the odd order derivatives are negative.
    I found this in the R-code of the equivalent filter.

MIGRATION TO python 3:
    Time to start migrating the codebase to python 3.

2017-03-31:
    In the fit_Hertz, Fmax can be None (default), thus no filtering will take 
    place. For now, the function will crash if F < Fmax is empty.
2017-04-04:
    Now, there is a general AFMforce class, which uses an internal dataset dict
    and extracts all information needed for the calculations.
    It has a dummy read_data() funciton, which is meant to be replaced to handle
    the input of specific AFM types.
    Now there is a JPKforce child and a Brukerforce child.
    The former handles the jpk-force binary files, the latter the text export
    from Bruker AFMs. I have tested with some Dimension Icon data so far.
    The test scripts are also renamed, to AFMIndentation, AFMPullAnalysis and
    AFMforcePreview. They can select which class they load either by config
    option 'datatype' or by extension:
    jpk-force or txt
    This is a great milestone of the project, making the system versatile. 
    Setup structure:
    The python files are now in the AFMforce folder, with proper __init__.py
    script, and in the Test folder as before.
    I plan further regrouping, so the classes sit in their files, and general
    analysis functions in their own set for better overview on the code.

2017-05-30:
    Baseline()... background correction. The current form shifts the i0 index
    to useless values, if the curve is well sampled. Thus, the local noise
    hides any trends, and the background is calculated from an adjusted,
    now wrong part of the curve.
    Change of strategy: do the fit to the user selected range, then calculate
    the error here. Extrapolate this error back, and find the maximum segment
    to the fit, which can be now refined.
    However, this refinement may be not necessary at all.

    Anyway, the code is modified to follow this strategy.

2017-10-06:
    A minor change. During fitting, we use only the x>d0 points for fitting
    however, this results in a flat line with F0 value if d0 > x.max() or
    d0 < x.min(). Thus, a flat area with moderate error that does not change
    at all. If this moderate error is small enough, the fitting ends with a
    warning that the error did not change, and that is it.
    Now I added a check to the fitting part for this cases to give some large
    value back hopefully pushing the fit back to track.

    In the AFMIndentation code I added the possibility to use the force curve
    from its minimum for the fit. This assumes that the attraction creates a 
    force offset only, but indentation happens from there already.
    The solution makes a pseudo curve keeping the data from the minimum, but
    padding it with some zeros to ensure a reasonable fit.

2018-05-11:
    Reworking the indentation part a bit to include surface tension terms,
    and esitmate contact behavior.

2018-05-14:
    An annoying bug in AFMIndentation: the Zmax parameter has to be applied
    after the selected segment was set! Setting the segment resets the Zmax
    property of the class, eliminating the previous settings!
    Thus, the correct way:
    - pick a file
    - select a segment
    - change Z-settings (e.g. Zmax, Z0, etc.)
    Actually, this is all wrong a bit.
    The invert() removes and reapplies Zmax. It should do it the correct way,
    then there should be no problem moving the segments around.
    Tested, it behaves the correct way now. I can swap between segments, and it
    does not change.

    O.k., this went a bit deeper than that.
    The point: because I used self.Zmax to decide if a filter was set, it broke
    down for cases where Zmax was set, but the segement was short.
    To clarify this issue, I switched checking if the _indx_ is empty,
    but leave Zmax alone. Why? Because intermediate segments of waiting at the
    surface or in retracted state have very short Z-ranges.
    With the current behavior, one sets Zmax once and it stays for all segments.

2018-06-25:
    Adding a gamma parameter to fit with gamma the Hertz model for spherical
    indenters. Here the point is to allow a surface effect, which is non-
    attractive. Now, sometimes the surface effect results in an unstable E-
    modulus. So, we can allow fixing the gamma (setting it non-negative), and
    fit only the Young's modulus and the surface zero point (zero indentation).
    This is now implemented, and the AFMIndentation test script can also handle
    it properly.
    It will result in E about 0 (-1E-5) if gamma is too large, because the data
    will become smaller than the lowest fitted curve. Naturally the fit would
    try going to negative E, which is penalized.

2018-07-24:
    I have trouble fititng. Many curves get an unreasonable high D0, and the
    whole gets shifted off....
    Solution: the automatic D0 estimation was way off. I did not activate this
    before, so it could not mess things up.
    Now it is way simpler than it was, and I think it is sufficiently so.
    The point: take a tail, and estimate standard deviation. Put Z0 to the
    lowest Z value where D < - std(D.end).
    And use the mean at the end as deflection 0.
    Why not the lowest Z-position of zero crossing? Because of the noise.
    Using -std() will be picking up the shift upwards earlier in noisy curves.

2018-12-07:
    I have trouble with SavGol, it does not preserve the peak position!
    Troubleshooting: use the R code as a reference to test what was going on.
    The nump.convolve is shifting the filtered points.
    Implement an fft based convolution, padding properly. The data is padded
    by its end value, the mask is by zeros as left to the fft itself.

2019-03-22:
    Still some problems. The JPK force data has noise both in Z and force
    directions. Now, if I smooth the Z data, I get jumps at the ends, because
    of the nonzero ends. I can remove one, but not both.
    Change the padding of the data before convolution to use the first and
    last points for the padding, and see if we get the curve back.
    The test SavGol filter works now.

    Now, I got the smoothing working right. The only thing to keep in mind,
    is that the kernels are not padded with their end values, or a boxcar
    kernel would be a huge constant. Thus, the provided kernels have to end
    with properly low values or they will shift the curve.
    For example, a Gaussian with 30,15 will make a shift, while with
    radius 30, width 5 will work like a charm.

    Next problem is handling JPK data:
    Here both Z and force are noisy. However, sorting Z evein in smoothened
    form will make the contact zone oscillating like hell.
    Use sorting only after identifying the contact zone and excluded it, or
    any peak finding will have million hits there.

    Now, we have the smoothing right. Because the noise in X and Y, this will
    shift the peaks a bit.
    I have tested the derivative filters. I have made Smooth capable of providing
    the first derivative for Savitzky - Golay and the Gauss filters. Using them
    to identify the minima has the problem, that their accuracy is related to
    the kernel width, and effects coming from the noise again.
    Best results I have got using a narraw, e.g. 5 pixel standard deviation
    Gaussian filter on the force data only. But even here the index identifycation
    is problematic, so at the end the running minimum filter is quite a nice
    solution it seems.


2019-11-21:
    A silly mistake: I always considered that the 'force' segment gets accessed
    with calibrated values. However, it is possible not to be the case.
    Thus, the part fishing for force constant gets called for converting raw
    values to voltages, where there is no distance data (deflection).
    In such a case an error occured.
    Now, that segment of code checks if there is a distance unit defined at all.

2020-01-29:
    Some fine tuning of the Gaussian filter. Working on fft based convolution,
    I learned that if the windo is too narrow, the nonzero offset will shift
    the fitted curve slightly down. It is more pronounced for curves with non-
    zero offset. And taking a note from the Numerical Recipes, I subtract the
    mean of the data before FFT. If I did smoothing (deriv==0), then I add it
    back after the convolution is complete. For calculating a derivative, this
    would shift the curve unnecessarily.
    I upped the verision to 0.62.

2021-01-27:
    Revisit the algorithms in working on a general analysis approach for force
    curves.
    Add smoothening parameters to the Contact_point function.

2021-02-11:
    Add an exponential fit. This was struggling until I realized that the x0
    parameter is useless here. It merges to the amplitude, thanks to the properties
    of the exponential.
    Now, it is working fine. The errors are super low though.

2021-03-09:
    Added the JKR fit. Here we are fitting the indentation to the F data, because
    the equation cannot be inverted. To improve the look, since the force data
    can be very noisy, an interpolated pair of variables as Fint and dint are also
    returned. These can be then plotted nicely to the curve.
    To plot the residuals, either sort the data first according to the x-axis used,
    or do not connect the points.

2021-03-26:
    Add position from the motorized stage to the piezo positions. This should allow
    for following other movements beyond the piexo stage, e.g. within the whole smple.

2021-05-19:
    Added extracting the date of the measurements from the files. This is a text, which
    may need conversion later on.

2021-05-25:
    added chi2 to the return of the Hertz model as well, since it saves recalculating it
    whenever needed. I have checked, it is exactly the same as doing by hand:
    ((res['y'] - res['fitted'])**2).sum()

2021-05-31:
    add a new fit function, which tries fitting a consecutive spherical Hertz model and
    an exponential to a force curve. It is best suited for cases where the 'contact' point
    is defined with the nonzero force (Contact_point use_deriv=False), and then see if a
    significant exponential contribution can be found in the end part.
    At the moment I cannot calculate the Jacoby matrix of the fitted function, because
    the index of the cut point is a free parameter in it. The fit converges almost to the
    suggested parameters. The first guess is based on a segmentwise linearized fit of the
    two functions.

2021-09-30:
    I am getting an error at line 1092 about min() being called on an empty array. I try
    getting around this if the indices are equal, use it as a result, but it may lead
    to other complications.

2022-06-20:
    Add a new class for Pavone nanoindenter measurements. This thing has a text export,
    with some special output. The sensor response is 1, because the data coming out is
    not really the same as for AFM.
    What we can do: we make the deflection set, raw data there is the same (this is why
    setting sensor response to 1, and calculate force from deflection (Cantilever).
    The measurement starts with Z=0 far away, so I turn the Z-axis around upon loading.

2022-06-24:
    Change to the setuptools using setup.cfg and a minimum pyproject.toml but no setup.py.
    This is not trivial, but using the documentation at:
    https://setuptools.pypa.io/en/latest/userguide/package_discovery.html
    Update:
        Using the setup.cfg from setuptools solved the problem. The pyproject.toml file
        just complicated the situation, but it did not work.

2023-03-01:
    Bugfix in the JPK reader. Sensor response may not be present, so define it only
    when we have it.

2023-08-25:
    The sloppy comparison of variable != [] for ndarrays does not work anymore, so I
    change the checks to variable.size > 0. Also make sure that whenever I have ndarrays,
    empty ones do not return [] --> lists, but nu.empty(0) an empty array.
    And remove the only left over nu.float... (deprecated)

2023-08-28:
    the slope detector (find_slope) has some issues time to time, so try fixing the code
    a bit. E.g. cleaning up the list comparisons, sorting, etc.

2024-04-20:
    Put the whole code under Creative Commons CC(4)-BY license, which allos free use,
    modifications, etc. but requires attribution.

2024-09-20:
    The Pavone has mixed up saving the retraction and contact segments, though they were
    registered correctly. I have added a force range test to try fixing it.

2025-05-02:
    Quite a bit happened recently. There are two key topics.
    1. an empyrical 7 parameter wave pattern fit is made to interference backgrounds
    2. the JPKforce part got some cleaning and update, because our new JPK system showed some
        quirks related to file format.

    The wavy background caused by light interference is not the simplest because it has:
    - periodicity that may shift linearly
    - amplitude that may increase / decrease
    - a tilt in the background
    So I made a triple 'linear' approximation of this problem adding amplitude slope,
    chirp and linear background to the curve.
    If the periodicity is quessed right enough it works. Otherwise it will settle to a
    zero amplitude middle value.
    This periodicity I guess from a smoothed data set and its crossing zeros.

    JPK data
    the company started storing the force/deflection/raw curves not in sequence but like
    raw/force/deflection, which meant sensor response was not picked up in time to get the
    force constant.
    But at the same time both of these calibration data are in a generic segment header
    copied in every segment.
    So, I added a sort of keys in evaluating the process and I am looking for the header as well.
    If I find it, it overwrites what I guessed from the slope values.

