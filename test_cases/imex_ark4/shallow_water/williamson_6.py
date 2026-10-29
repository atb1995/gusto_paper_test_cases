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
    implicit, explicit, split_continuity_form, TrapeziumRule, BackwardEuler,
    IMEX_ARS3, IMEX_SSP3, IMEX_Euler, IMEXRungeKutta
)
import time
import numpy as np

williamson_6_defaults = {
    'ncells_per_edge': 32,     # number of cells per cubed sphere panel edge
    'dt': 900.0,               # 15 minutes
    'tmax': 14.*24.*60.*60.,   # 14 days
    'dumpfreq': 96,            # once per day with default options
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

    # Make sure divergence term is split
    eqns = split_continuity_form(eqns)
    transport_methods = [
        DGUpwind(eqns, "u"),
        DGUpwind(eqns, "D")
    ]

    # Hybridised solver settings
    settings = {
            'snes_type': 'newtonls',
            'snes_max_it': 50,
            'snes_atol': 1e-4,
            'snes_rtol': 1e-4,
            
            'snes_monitor': None,
            'snes_converged_reason': None,
            'ksp_type': 'preonly',
            'mat_type': 'matfree',
            'pc_type': 'python',
            'pc_python_type': 'firedrake.HybridizationPC',  # Uses Firedrake's
                'hybridization': {                              # hybridization PC
                    'ksp_type': 'gmres',
                    'ksp_max_it': 200,
                    'pc_type': 'gamg',  # AMG for trace system
                    'ksp_rtol': 1e-8,
                    'mg_levels': {
                        'ksp_type': 'chebyshev',
                        'ksp_max_it': 2,
                        'pc_type': 'bjacobi',
                        'sub_pc_type': 'ilu'
                    }
                }
        }
    
    # Transport terms explicit, all other terms implicit
    eqns.label_terms(lambda t: t.has_label(transport), explicit)
    eqns.label_terms(lambda t:  not t.has_label(time_derivative) and not t.has_label(explicit), implicit)
    
    # IMEX-ARK4 butcher tables
    butcher_imp=np.array(
        [
            [0, 0, 0, 0, 0, 0],
            [1 / 4, 1 / 4, 0, 0, 0, 0],
            [8611 / 62500, -1743 / 31250, 1 / 4, 0, 0, 0],
            [5012029 / 34652500, -654441 / 2922500, 174375 / 388108, 1 / 4, 0, 0],
            [
                15267082809 / 155376265600,
                -71443401 / 120774400,
                730878875 / 902184768,
                2285395 / 8070912,
                1 / 4,
                0,
            ],
            [82889 / 524892, 0, 15625 / 83664, 69875 / 102672, -2260 / 8211, 1 / 4],
            [82889 / 524892, 0, 15625 / 83664, 69875 / 102672, -2260 / 8211, 1 / 4]
        ],
        dtype=float,
    )
    butcher_exp=np.array(
        [
            [0, 0, 0, 0, 0, 0],
            [1 / 2, 0, 0, 0, 0, 0],
            [13861 / 62500, 6889 / 62500, 0, 0, 0, 0],
            [
                -116923316275 / 2393684061468,
                -2731218467317 / 15368042101831,
                9408046702089 / 11113171139209,
                0,
                0,
                0,
            ],
            [
                -451086348788 / 2902428689909,
                -2682348792572 / 7519795681897,
                12662868775082 / 11960479115383,
                3355817975965 / 11060851509271,
                0,
                0,
            ],
            [
                647845179188 / 3216320057751,
                73281519250 / 8382639484533,
                552539513391 / 3454668386233,
                3354512671639 / 8306763924573,
                4040 / 17871,
                0,
            ],
            [82889 / 524892, 0, 15625 / 83664, 69875 / 102672, -2260 / 8211, 1 / 4]
            
        ],
        dtype=float,
    )

    # IMEX-ARK4 timestepper
    scheme = IMEXRungeKutta(domain, butcher_imp=butcher_imp, butcher_exp=butcher_exp, nonlinear_solver_parameters=settings)
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
    scheme_name = "IMEX_ARK5"
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
