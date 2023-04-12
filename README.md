# AFM force package for python
A set of functions to load force-distance data from JPK and Bruker AFMs, the latter
only from ASCII exported files.

## Aim
load and analyze force-distance curves in a relatively easy way.

## Main class: AFMforce
An envelop to contain piezo actuator distance data and deflection data detected
in the instrument. Internal calibration constants allow then to turn the sensor
data to deflection and force.
The main class has variants for the various instruments, providing a unified
interface to the user for studying the curves.

# Units
Units are provided to guide the user what is provided. E.g. JPK usually uses
direct SI units: m for distance and N for force. Others may use nm and nN
respectively.
The analysis functions assume distances in micrometers, forces in nN.

# Smoothing
A Savitzky - Golay filter or a Gaussian smoothing is available when
needed. The non-zero padding of the convolution remvoes the cut-off
effect of common convolution, making both filters working excellent on the data.

# fitting
## Indentation
The standard Hertz model is implemented for conical tips or spherical
probes. The algorithm is a simple nonlinear fit using the leastsq()
algorithm from scipy. Some linearized guessing is done to have a reasonable
starting point of the routine.

### gamma
An alternative fit is provided where a line tension term is also included,
adding a straight line trend to the fit. The gamma part can be either fixed
or a fit parameter itself. This latter is often results in a weaker fit
behavior, so be careful when using it.

## exponential
Force distance curves can be also fit to an exponential decay, which
may be the result of chain overlapping of polymer brushes of a DLVO-like
interaction between surfaces. A zero point is meaningless for an exponential
curve, because it automatically incorporates into the amplitude of the
curve.

## Hertz - Exp
A segmented fit is also available, which assumes an exponential decay
at higher distances, then a Hertzian curve. It tries matching the two
segments, but setting the crossing point a free parameter turns the
process to an ill-conditioned one.

## JKR
for the case of adhesion between the surface and the probe, a JKR model
fit is also provided.

## Pull
Analysing pull away events, a worm-like chain model based fit is available.
Though the fit works, the resulted parameters are rather sensitive to the
experimental conditions.

## LM
a generic linear model is available to combine any data to a linear
combination of known functions. It is a similar algorithm to that of
the lm() function in R. (And was written before numpy got anything
similar on the line 8).)

# base line and force vs. distance
A base line fit is provided to estimate the zero force background
and to help obtaining force - distance curves.
Another algorithm is available to find a contact point.
If there are attraction forces, then the contact point is
where the retract curve crosses the zero force base line.
If not, one can dfine a contact point where the force deviates from
zero base line more than 3x the noise level (standard deviation of
the zero force base line).

(Both assumptions may be acceptable or may fail time to time.)

# Other functions
There are functions to estimate whether puncture events happened
during an approach, or to find what is the hard wall slope of
a curve.

FindMinima is to find local minima identifying break-away events
in pull-back curves.
