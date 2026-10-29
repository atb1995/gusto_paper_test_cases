"""
The Rossby-Haurwitz Wave test case (6) of Williamson et al. 1992.
A non-divergent wave pattern moves eastwardly, maintaining the wave structure.
The test case is on the sphere.

The setup implemented here uses the cubed sphere mesh with the degree 1 spaces.
"""
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter

from firedrake import (
    SpatialCoordinate, Constant,
    Function, sin, cos, FunctionSpace, grad
)
from gusto import (
    Domain, IO, OutputParameters, GeneralCubedSphereMesh, RelativeVorticity,
    lonlatr_from_xyz, ShallowWaterEquations, ShallowWaterParameters, SSPRK3,
    DGUpwind, SemiImplicitQuasiNewton, ZonalComponent, MeridionalComponent,
    SubcyclingOptions, RungeKuttaFormulation, Timestepper, divergence, 
    pressure_gradient, gravity, coriolis, transport, time_derivative,
    implicit, explicit, split_continuity_form, TrapeziumRule, BackwardEuler
)
import time

williamson_6_defaults = {
    'ncells_per_edge': 128,     # number of cells per cubed sphere panel edge
    'dt': 3.0,               # 15 minutes
    'tmax': 1.*24.*60.*60.,   # 1 day
    'dumpfreq': 14400,            # once per day with default options
    'dirname': 'williamson_6'
}


def williamson_6(
        ncells_per_edge=williamson_6_defaults['ncells_per_edge'],
        dt=williamson_6_defaults['dt'],
        tmax=williamson_6_defaults['tmax'],
        dumpfreq=williamson_6_defaults['dumpfreq'],
        dirname=williamson_6_defaults['dirname']
):

    # ------------------------------------------------------------------------ #
    # Test case parameters
    # ------------------------------------------------------------------------ #

    radius = 6371220.       # radius of the planet, in m
    K = Constant(7.847e-6)  # Frequency parameter, in s^{-1}
    w = K                   # Set omega equal to K
    R = 4.                  # Wave number 4
    h0 = 8000.              # mean (and reference) depth, in m

    # ------------------------------------------------------------------------ #
    # Our settings for this set up
    # ------------------------------------------------------------------------ #

    degree = 1
    hdiv_family = 'RTCF'
    u_eqn_type = 'vector_advection_form'

    # ------------------------------------------------------------------------ #
    # Set up model objects
    # ------------------------------------------------------------------------ #

    # Domain
    mesh = GeneralCubedSphereMesh(radius, ncells_per_edge, degree=2)
    domain = Domain(mesh, dt, hdiv_family, degree)

    # Equation
    xyz = SpatialCoordinate(mesh)
    parameters = ShallowWaterParameters(mesh, H=h0)
    eqns = ShallowWaterEquations(
        domain, parameters, u_transport_option=u_eqn_type
    )

    # I/O and diagnostics
    output = OutputParameters(
        dirname=dirname, dumpfreq=dumpfreq, dump_nc=True, dumplist=['u', 'D'],
        checkpoint=True, checkpoint_method="checkpointfile", chkptfreq=dumpfreq
    )
    diagnostic_fields = [
        RelativeVorticity(), ZonalComponent('u'), MeridionalComponent('u')
    ]
    io = IO(domain, output, diagnostic_fields=diagnostic_fields)

    # Transport schemes
    transport_methods = [
        DGUpwind(eqns, "u"),
        DGUpwind(eqns, "D")
    ]
    # SSPRK3 Timestepper
    scheme = SSPRK3(domain)
    stepper = Timestepper(eqns, scheme, io, transport_methods )

    # ------------------------------------------------------------------------ #
    # Initial conditions
    # ------------------------------------------------------------------------ #

    u0_field = stepper.fields("u")
    D0_field = stepper.fields("D")

    # Parameters
    g = parameters.g
    Omega = parameters.Omega

    lon, lat, _ = lonlatr_from_xyz(xyz[0], xyz[1], xyz[2])

    # Initialise the velocity field from streamfunction
    CG2 = FunctionSpace(mesh, 'CG', 2)
    psi = Function(CG2)
    psiexpr = (
        -radius**2 * w * sin(lat)
        + radius**2 * K * cos(lat)**R * sin(lat) * cos(R*lon)
    )
    psi.interpolate(psiexpr)
    uexpr = domain.perp(grad(psi))

    # Depth field
    A = (
        (w / 2) * (2 * Omega + w) * cos(lat)**2
        + 0.25 * K**2 * cos(lat)**(2 * R) * (
            (R + 1) * cos(lat)**2 + (2*R**2 - R - 2) - 2*R**2 * cos(lat)**(-2)
        )
    )
    B_frac = (2 * (Omega + w) * K) / ((R + 1) * (R + 2))
    B = B_frac * cos(lat)**R * ((R**2 + 2 * R + 2) - (R + 1)**2 * cos(lat)**2)
    C = (1 / 4) * K**2 * cos(lat)**(2 * R) * ((R + 1)*cos(lat)**2 - (R + 2))
    Dexpr = h0 * g + radius**2 * (A + B*cos(lon*R) + C * cos(2 * R * lon))

    # Obtain fields
    u0_field.project(uexpr)
    D0_field.interpolate(Dexpr / g)

    # Set reference profiles
    Dbar = Function(D0_field.function_space()).assign(D0_field)
    stepper.set_reference_profiles([('D', Dbar)])

    # ------------------------------------------------------------------------ #
    # Run
    # ------------------------------------------------------------------------ #
    start_time = time.time()
    stepper.run(t=0, tmax=tmax)
    end_time = time.time()
    scheme_name = "SSPRK3"
    dt_val = dt
    tmax_val = tmax
    ncells = ncells_per_edge
    elapsed_s = end_time - start_time
    print(
        "RESULT,"
        f"scheme={scheme_name},"
        f"dt_s={dt_val:.6f},"
        f"tmax_s={tmax_val:.6f},"
        f"ncells_per_edge={ncells},"
        f"elapsed_s={elapsed_s:.6f},"
    )


# ---------------------------------------------------------------------------- #
# MAIN
# ---------------------------------------------------------------------------- #


if __name__ == "__main__":

    parser = ArgumentParser(
        description=__doc__,
        formatter_class=ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        '--ncells_per_edge',
        help="The number of cells per edge of cubed sphere panel",
        type=int,
        default=williamson_6_defaults['ncells_per_edge']
    )
    parser.add_argument(
        '--dt',
        help="The time step in seconds.",
        type=float,
        default=williamson_6_defaults['dt']
    )
    parser.add_argument(
        "--tmax",
        help="The end time for the simulation in seconds.",
        type=float,
        default=williamson_6_defaults['tmax']
    )
    parser.add_argument(
        '--dumpfreq',
        help="The frequency at which to dump field output.",
        type=int,
        default=williamson_6_defaults['dumpfreq']
    )
    parser.add_argument(
        '--dirname',
        help="The name of the directory to write to.",
        type=str,
        default=williamson_6_defaults['dirname']
    )
    args, unknown = parser.parse_known_args()

    williamson_6(**vars(args))
