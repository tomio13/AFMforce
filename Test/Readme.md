# Test
These are not so much classical test cases, more like
try using the library to do some interesting automated
analysis of force curves.

# AFMforcePreview
a simple script to read any known force formats, and provide
a simple matplotlib plot on screen. It helps identifying really
ill formed curves to be trashed before running a batch process.


# AFMIndentation
run a Hertz model on the curves. It can handle approach or
retraction sides, and a couple of parameters.

Many of these scripts use the BatchAnalyzer library (of like
5 functions) to simpler configuration and reporting.
It is a super simple thing to do some easy handling, and will
be uploaded to Github for usage. But do not expect anything
fancy 8).

# AFM-force-analysis
run a generic analysis trying to find exponential and spherical
Hertzian indentation segments, see if the generic indentaion
matches the data, find maximal adhesion force, adhesion energy,
contact thickness, etc. A bit of heuristic treatment, but it
is often useful to have a generic idea of the force distance
curves.
