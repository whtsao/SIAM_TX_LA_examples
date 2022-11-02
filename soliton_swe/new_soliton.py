from __future__ import division
from builtins import object
from past.utils import old_div
from proteus.mprans import (SW2DCV, GN_SW2DCV)
from proteus.Domain import RectangularDomain, PlanarStraightLineGraphDomain
import numpy as np
from proteus import (Domain, Context, MeshTools as mt)
from proteus.Profiling import logEvent
from proteus.Gauges import PointGauges
import proteus.SWFlow.SWFlowProblem as SWFlowProblem
import os

"""
This is the problem of a solitary wave run-up on a sloping beach.
"""

# *************************** #
# ***** GENERAL OPTIONS ***** #
# *************************** #
opts = Context.Options([
    # Simulation options
    ('sw_model', 1, "sw_model = {0,1} for {SWEs,DSWEs}"),
    ("final_time", 60., "Final time for simulation"),
    ("dt_output", 0.1, "Time interval to output solution"),
    ("cfl", 0.25, "Desired CFL restriction"),
    # Mesh options
    ("refinement", 6, "Refinement level"),
    ("structured", True, "Structured or unstructured mesh"),
    ("he", 0.5, "Mesh size for unstructured mesh"),
    ("reflecting_BCs", False, "Use reflecting BCs for all boundaries"),
    # Problem specific options
    ("want_gauges", True, "Output for water height point gauge"),
    ("mannings", 0., "Mannings roughness coefficient"), # usually = 0
    ("still_water_depth", 0.7, "Depth of still water above floor"),
    ("solitary_amplitude", 0.07, "Amplitude of solitary wave"),
    ("solitary_position", 10., "Center position of circular dam"),
    ("ramp_slope", 1. / 50., "Slope of ramp"),
    ("ramp_location", 55., "Location where ramp starts")
])

###################
# DOMAIN AND MESH #
###################
L = (75., 4.)  # this is length in x direction and y direction
refinement = opts.refinement
origin = [0, 0, 0]
rectangle = RectangularDomain(L=L, x=origin)

# CREATE REFINEMENT #
nnx0 = 6
nnx = (nnx0 - 1) * (2**refinement) + 1
nny = old_div((nnx - 1), 10) + 1
he = old_div(L[0], float(nnx - 1))
if opts.structured:
    domain = rectangle
else:
    rectangle.writePoly("sloping_beach")
    domain = PlanarStraightLineGraphDomain(fileprefix="sloping_beach")
    domain.MeshOptions.triangleOptions = "pAq30Dena%f" % (0.5 * opts.he**2,)
    nnx = None
    nny = None

###############################
#  CONSTANTS NEEDED FOR SETUP #
###############################
g = 9.81

# constants solitary wave
h0 = opts.still_water_depth
alpha = opts.solitary_amplitude
xs = opts.solitary_position
r = np.sqrt(old_div(3. * alpha, (4. * h0**2 * (h0 + alpha))))
c = np.sqrt(g * (h0 + alpha))

# constants for bathymetry
ramp_slope = opts.ramp_slope
ramp_location = opts.ramp_location


#####################################
#   Some functions defined here    #
####################################

def solitary_wave(x, t):
    sechSqd = (1.0 / np.cosh(r * (x - xs - c * t)))**2
    return alpha * sechSqd


def bathymetry_function(X):
    x = X[0]
    y = X[1]

    # first define sloping ramp
    ramp = ramp_slope * (x - ramp_location)

    # define base
    base = 0. * x

    # define final bathymetry
    bath = np.maximum(base, ramp)
    return bath


######################
# INITIAL CONDITIONS #
######################


class water_height_at_t0(object):
    def uOfXT(self, X, t):
        hTilde = h0 + solitary_wave(X[0], 0)
        h = max(hTilde - bathymetry_function(X), 0.)
        return h


class x_mom_at_t0(object):
    def uOfXT(self, X, t):
        hTilde = h0 + solitary_wave(X[0], 0)
        h = max(hTilde - bathymetry_function(X), 0.)
        return h * c * old_div(hTilde - h0, hTilde)

"""
heta and hw are needed for the hyperbolic serre-green-naghdi equations.
For initial conditions, heta -> h^2, hbeta->q(dot)grad(Z), hw -> h^2div(u)+3/2*hbeta.
It's often okay to take hbeta=0. Note that the BCs for the heta and hw should be same as h
and BCs for hbeta should be same as x_mom.
For more details see: 'Hyperbolic relaxation technique for solving the dispersive Serre--Green-Naghdi Equations
with topography' by Guermond, Kees, Popov, Tovar.
"""

class heta_at_t0(object):
    def uOfXT(self, X, t):
        h = water_height_at_t0().uOfXT(X, t)
        return h**2


class hw_at_t0(object):
    def uOfXT(self, X, t):
        sechSqd = (1.0 / np.cosh(r * (X[0] - xs)))**2.0
        hTilde = h0 + solitary_wave(X[0], 0)
        h = max(hTilde - bathymetry_function(X), 0.)
        hTildePrime = -2.0 * alpha * r * np.tanh(r * (X[0] - xs)) * sechSqd
        hw = -h**2 * old_div(c * h0 * hTildePrime, hTilde**2)
        return hw

class Zero(object):
    def uOfXT(self, X, t):
        return 0.

###############################
##### BOUNDARY CONDITIONS #####
###############################

# this is x domain, used in BCs
X_coords = (origin[0], origin[0] + L[0])

# this is y domain, used in BCs
Y_coords = (origin[1], origin[1] + L[1])


def x_mom_DBC(X, flag):
    if X[0] == X_coords[0] or X[0] == X_coords[1]:
        return lambda X, t: 0.0


def y_mom_DBC(X, flag):
    if X[1] == Y_coords[0] or X[1] == Y_coords[1]:
        return lambda X, t: 0.0


# ********************************** #
# ***** Create mySWFlowProblem ***** #
# ********************************** #
outputStepping = SWFlowProblem.OutputStepping(
    opts.final_time, dt_output=opts.dt_output)
initialConditions = {'water_height': water_height_at_t0(),
                     'x_mom': x_mom_at_t0(),
                     'y_mom': Zero(),
                     'h_times_eta': heta_at_t0(),
                     'h_times_w': hw_at_t0(),
                     'h_times_beta': Zero()}
boundaryConditions = {'water_height': lambda x, flag: None,
                      'x_mom': x_mom_DBC,
                      'y_mom': y_mom_DBC,
                      'h_times_eta': lambda x, flag: None,
                      'h_times_w': lambda x, flag: None,
                      'h_times_beta': x_mom_DBC}
# **************************** #
# ********** GAUGES ********** #
# **************************** #

# x location of gauges
first_location = 55.0
second_location = first_location + 16.25
third_location = second_location + 1.5

# y location of gauges
mid_of_domain = 0.5 * L[1]


heightPointGauges = PointGauges(gauges=((('h',),
                                ((first_location, mid_of_domain, 0), (second_location, mid_of_domain, 0),
                                (third_location, mid_of_domain, 0))),),
                                activeTime=(0., opts.final_time),
                                fileName='beach_wave_gauges.csv')

# ********************************************* #
# ********** Create my SWFlowProblem ********** #
# ********************************************* #
mySWFlowProblem = SWFlowProblem.SWFlowProblem(sw_model=opts.sw_model,
                                              cfl=opts.cfl,
                                              outputStepping=outputStepping,
                                              structured=opts.structured,
                                              he=he,
                                              nnx=nnx,
                                              nny=nny,
                                              domain=domain,
                                              initialConditions=initialConditions,
                                              boundaryConditions=boundaryConditions,
                                              reflectingBCs=opts.reflecting_BCs,
                                              bathymetry=bathymetry_function,
                                              analyticalSolution=None)
mySWFlowProblem.physical_parameters['mannings'] = opts.mannings
if opts.want_gauges:
    mySWFlowProblem.auxiliaryVariables = [heightPointGauges]
