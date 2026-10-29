"""
The deep atmosphere dry baroclinic wave from Ullrich et al. 2014:
``A proposed baroclinic wave test case for deep- and shallow-atmosphere
dynamical cores'', QJRMS.

This is a 3D test on the sphere, with an initial state that is in unsteady
balance, with a perturbation added to the wind.

This setup uses a cubed-sphere with the order 1 finite elements.
"""
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter, ArgumentTypeError

from firedrake import (
    ExtrudedMesh, SpatialCoordinate, cos, sin, pi, sqrt, exp, Constant,
    Function, acos, errornorm, norm, le, ge, conditional, inner, dx,
    NonlinearVariationalProblem, NonlinearVariationalSolver, TestFunction,
    min_value, max_value
)
from gusto import *
from mpi4py import MPI

import csv
import os
import time

# ---------------------------------------------------------------------------- #
# Test case defaults
# ---------------------------------------------------------------------------- #

def create_layer_height(htop, nlayers):
    """
    Compute the layer heights for the dry baroclinic sphere test case.

    Parameters
    ----------
    htop : float
        The height of the top of the atmosphere above the surface, in m.
    nlayers : int
        The number of layers in the mesh.

    Returns
    -------
    list of float
        The heights of each layer in the mesh.
    """
    layer_height = []
    running_height = 0
    for m in range(1, nlayers+1):
        mu = 15
        height = htop * (
            ((mu * (m / nlayers)**2 + 1)**0.5 - 1)
            / ((mu + 1)**0.5 - 1)
        )
        depth = height - running_height
        running_height = height
        layer_height.append(depth)

    return layer_height

dry_baroclinic_sphere_defaults = {
    'ncell_per_edge': 16,
    'nlayers': 15,
    'dt': 450.0,               # 15 minutes
    'tmax': 15.*24.*60.*60.,   # 15 days
    'dumpfreq': 100,            # Corresponds to every 12 hours with default opts
    'dirname': 'dry_baroclinic_sphere_siqn',
    'results_csv': None,
    'io_on':True
}


def dry_baroclinic_sphere(
        ncell_per_edge=dry_baroclinic_sphere_defaults['ncell_per_edge'],
        nlayers=dry_baroclinic_sphere_defaults['nlayers'],
        dt=dry_baroclinic_sphere_defaults['dt'],
        tmax=dry_baroclinic_sphere_defaults['tmax'],
        dumpfreq=dry_baroclinic_sphere_defaults['dumpfreq'],
        dirname=dry_baroclinic_sphere_defaults['dirname'],
        results_csv=dry_baroclinic_sphere_defaults['results_csv'],
        io_on=dry_baroclinic_sphere_defaults['io_on']
):

    # ------------------------------------------------------------------------ #
    # Parameters for test case
    # ------------------------------------------------------------------------ #
    a = 6.371229e6    # radius of planet, in m
    htop = 3.0e4      # height of top of atmosphere above surface, in m
    omega = 7.292e-5  # rotation frequency of planet, in 1/s
    lapse = 0.005     # thermal lapse rate
    T0e = 310         # equatorial temperature, in K
    T0p = 240         # polar surface temperature, in K
    Vp = 1.0          # maximum velocity of perturbation in m/s
    z_pert = 1.5e4    # height of perturbation
    d0 = a/6          # horizontal radius of perturbation
    lon_c = pi/9      # longitude of perturbation centre
    lat_c = 2*pi/9    # latitude of perturbation centre

    # ------------------------------------------------------------------------ #
    # Our settings for this set up
    # ------------------------------------------------------------------------ #
    horder = 1        # horizontal order of finite element de Rham complex
    vorder = 1        # vertical order of finite element de Rham complex
    u_eqn_type = 'vector_advection_form'  # Form of the momentum equation to use

    # ------------------------------------------------------------------------ #
    # Set up model objects
    # ------------------------------------------------------------------------ #

    # Layers are not evenly spaced
    layer_height = create_layer_height(htop, nlayers)

    from firedrake import DistributedMeshOverlapType

    distribution_parameters = {
        "partition": True,
        "overlap_type": (DistributedMeshOverlapType.VERTEX, 1)
    }

    base_mesh = GeneralCubedSphereMesh(
        a, num_cells_per_edge_of_panel=ncell_per_edge, degree=2,
        distribution_parameters=distribution_parameters
    )
    mesh = ExtrudedMesh(
        base_mesh, layers=nlayers, layer_height=layer_height,
        extrusion_type='radial'
    )
    domain = Domain(
        mesh, dt, "RTCF", horizontal_degree=horder, vertical_degree=vorder
    )
    params = CompressibleParameters(mesh, Omega=omega)
    eqn = CompressibleEulerEquations(
        domain, params, u_transport_option=u_eqn_type
    )

    print("Opt Cores:", eqn.X.function_space().dim()/50000.)

    # ------------------------------------------------------------------------ #
    # IO
    # ------------------------------------------------------------------------ #
    if not io_on:
        dumpfreq = 999
        dump_nc = False
        dump_vtus = False
        dumplist = []

        diagnostic_fields = []
    else:
        dump_nc = True
        dump_vtus = False
        dumplist = ['rho', 'theta']
        diagnostic_fields = [
            ZonalComponent('u'), MeridionalComponent('u'),
            Temperature(eqn), Pressure(eqn)
        ]
    output = OutputParameters(
        dirname=dirname, dumpfreq=dumpfreq, dump_nc=dump_nc, dump_vtus=dump_vtus, dumplist=dumplist
    )
    io = IO(domain, output, diagnostic_fields=diagnostic_fields)

    # ------------------------------------------------------------------------ #
    # Transport
    # ------------------------------------------------------------------------ #

    theta_opts = EmbeddedDGOptions()
    transported_fields = [
        SSPRK3(domain, "u"),
        SSPRK3(domain, "rho", rk_formulation=RungeKuttaFormulation.linear),
        SSPRK3(domain, "theta", options=theta_opts)
    ]
    transport_methods = [
        DGUpwind(eqn, "u"),
        DGUpwind(eqn, "rho", advective_then_flux=True),
        DGUpwind(eqn, "theta")
    ]

    # ------------------------------------------------------------------------ #
    # Time stepper
    # ------------------------------------------------------------------------ #

    stepper = SemiImplicitQuasiNewton(
        eqn, io, transported_fields, transport_methods, predictor='rho',
        tau_values={'rho': 1.0, 'theta': 1.0}
    )

    # ------------------------------------------------------------------------ #
    # Initial conditions
    # ------------------------------------------------------------------------ #

    x, y, z = SpatialCoordinate(mesh)
    lon, lat, r = lonlatr_from_xyz(x, y, z)

    u0 = stepper.fields("u")
    rho0 = stepper.fields("rho")
    theta0 = stepper.fields("theta")

    Vu = u0.function_space()
    Vr = rho0.function_space()
    Vt = theta0.function_space()

    # Steady state -------------------------------------------------------------

    Rd = params.R_d
    g = params.g
    p0 = params.p_0

    T0 = 0.5*(T0e + T0p)
    H = Rd*T0/g
    k = 3
    b = 2

    s = (r/a)*cos(lat)
    A = 1/lapse
    B = (T0e - T0p)/((T0e + T0p)*T0p)
    C = ((k + 2)/2)*((T0e - T0p)/(T0e*T0p))

    tau1 = A*lapse*exp((r - a)*lapse/T0)/T0
    tau1 += B * (1 - 2*((r - a)/(b*H))**2) * exp(-((r - a) / (b*H))**2)

    tau2 = C * (1 - 2*((r - a)/(b*H))**2) * exp(-((r - a) / (b*H))**2)

    tau1_integral = A * (exp(lapse * (r - a) / T0) - 1)
    tau1_integral += B * (r - a) * exp(-((r - a) / (b*H))**2)

    tau2_integral = C * (r - a) * exp(-((r - a) / (b*H))**2)

    T_expr = (a / r)**2 / (
        tau1 - tau2 * (s**k - (k/(k + 2)) * s**(k + 2))
    )
    P_expr = p0 * exp(
        - g/Rd * tau1_integral
        + g/Rd * tau2_integral * (s**k - (k / (k + 2)) * s**(k + 2))
    )

    wind_proxy = (
        (g/a)*k*T_expr*tau2_integral*(
            ((r*cos(lat))/a)**(k - 1) - ((r*cos(lat))/a)**(k + 1)
        )
    )
    wind = (
        - omega*r*cos(lat)
        + sqrt((omega*r*cos(lat))**2 + r*cos(lat)*wind_proxy)
    )

    theta_expr = T_expr*(P_expr/p0)**(- params.kappa)
    exner_expr = T_expr/theta_expr
    rho_expr = P_expr/(Rd*T_expr)

    # Perturbation -------------------------------------------------------------

    acos_arg = sin(lat_c)*sin(lat) + cos(lat_c)*cos(lat)*cos(lon - lon_c)
    acos_arg_clamped = min_value(1.0, max_value(-1.0, acos_arg))
    d = a * acos(acos_arg_clamped)

    height = r - a
    err_tol = 1e-12
    zeta = conditional(
        ge(height, z_pert - err_tol),
        0,
        1 - 3*(height/z_pert)**2 + 2*(height/z_pert)**3
    )

    perturb_magnitude = (
        (16*Vp/(3*sqrt(3)))*zeta*sin((pi*d)/(2*d0))*cos((pi*d)/(2*d0))**3
    )

    zonal_pert = conditional(
        le(d, err_tol),
        0,
        conditional(
            ge(d, (d0 - err_tol)),
            0,
            - perturb_magnitude*(
                -sin(lat_c)*cos(lat) + cos(lat_c)*sin(lat)*cos(lon - lon_c)
            )/sin(d/a)
        )
    )
    meridional_pert = conditional(
        le(d, err_tol),
        0,
        conditional(
            ge(d, d0 - err_tol),
            0,
            perturb_magnitude*cos(lat_c)*sin(lon - lon_c)/sin(d/a)
        )
    )

    zonal_u = wind + zonal_pert
    merid_u = Constant(0.0) + meridional_pert
    radial_u = Constant(0.0)

    e_lon = xyz_vector_from_lonlatr(1, 0, 0, (x, y, z))
    e_lat = xyz_vector_from_lonlatr(0, 1, 0, (x, y, z))
    e_r = xyz_vector_from_lonlatr(0, 0, 1, (x, y, z))

    test_u = TestFunction(Vu)
    dx_reduced = dx(degree=8)
    u_field = zonal_u*e_lon + merid_u*e_lat + radial_u*e_r
    u_proj_eqn = inner(test_u, u0)*dx - inner(test_u, u_field)*dx_reduced
    u_proj_prob = NonlinearVariationalProblem(u_proj_eqn, u0, bcs=eqn.bcs['u'])
    u_solver_parameters = {
        'snes_type': 'ksponly',
        'ksp_rtol': 1e-12,
        'ksp_type': 'cg',
        'ksp_max_it': 200,
        'ksp_monitor_true_residual': None,
        'ksp_converged_reason': None,
        'pc_type': 'hypre',
    }
    u_proj_solver = NonlinearVariationalSolver(u_proj_prob, solver_parameters=u_solver_parameters)
    u_proj_solver.solve()

    theta0.interpolate(theta_expr)

    exner = Function(Vr).interpolate(exner_expr)
    rho0.interpolate(rho_expr)

    compressible_hydrostatic_balance(
        eqn, theta0, rho0, exner_boundary=exner, solve_for_rho=True
    )

    rho_analytic = Function(Vr).interpolate(rho_expr)
    print(f'Normalised rho error: {errornorm(rho_analytic, rho0)/norm(rho_analytic)}')

    rho_b = Function(Vr).assign(rho0)
    theta_b = Function(Vt).assign(theta0)
    stepper.set_reference_profiles([('rho', rho_b), ('theta', theta_b)])

    # ------------------------------------------------------------------------ #
    # Run
    # ------------------------------------------------------------------------ #

    start_time = MPI.Wtime()
    stepper.run(t=0, tmax=tmax)
    end_time = MPI.Wtime()
    first_step_time = stepper.first_step_time
    wall_time = end_time - start_time
    first_step_wall_time = first_step_time - start_time
    from_first_step_wall_time = wall_time - first_step_wall_time
    print(f"First step complete. dt={dt}, wall time={first_step_wall_time:.1f}s")
    print(f"Run complete. dt={dt}, wall time={wall_time:.1f}s")

    comm = MPI.COMM_WORLD
    if comm.Get_rank() == 0 and results_csv is not None:
        write_header = not os.path.exists(results_csv)
        os.makedirs(os.path.dirname(results_csv) or ".", exist_ok=True)
        with open(results_csv, 'a', newline='') as f:
            writer = csv.writer(f)
            if write_header:
                writer.writerow(['ncells', 'nlayers', 'dt', 'tmax', 'nprocs', 'first_step_wall_time', 'wall_time', 'from_first_step_wall_time'])
            writer.writerow([ncell_per_edge, nlayers, dt, tmax, comm.Get_size(), f"{first_step_wall_time:.1f}", f"{wall_time:.1f}", f"{from_first_step_wall_time:.1f}"])


def parse_bool(value):
    if isinstance(value, bool):
        return value
    value = str(value).strip().lower()
    if value in {"1", "true", "t", "yes", "y", "on"}:
        return True
    if value in {"0", "false", "f", "no", "n", "off"}:
        return False
    raise ArgumentTypeError(f"Boolean value expected, got {value!r}")


# ---------------------------------------------------------------------------- #
# MAIN
# ---------------------------------------------------------------------------- #

if __name__ == "__main__":

    parser = ArgumentParser(
        description=__doc__,
        formatter_class=ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        '--ncell_per_edge',
        help="The number of cells per panel edge of the cubed-sphere.",
        type=int,
        default=dry_baroclinic_sphere_defaults['ncell_per_edge']
    )
    parser.add_argument(
        '--nlayers',
        help="The number of layers for the mesh.",
        type=int,
        default=dry_baroclinic_sphere_defaults['nlayers']
    )
    parser.add_argument(
        '--dt',
        help="The time step in seconds.",
        type=float,
        default=dry_baroclinic_sphere_defaults['dt']
    )
    parser.add_argument(
        "--tmax",
        help="The end time for the simulation in seconds.",
        type=float,
        default=dry_baroclinic_sphere_defaults['tmax']
    )
    parser.add_argument(
        '--dumpfreq',
        help="The frequency at which to dump field output.",
        type=int,
        default=dry_baroclinic_sphere_defaults['dumpfreq']
    )
    parser.add_argument(
        '--dirname',
        help="The name of the directory to write to.",
        type=str,
        default=dry_baroclinic_sphere_defaults['dirname']
    )
    parser.add_argument(
        '--results_csv',
        help="Path to CSV file to append timing results to.",
        type=str,
        default=dry_baroclinic_sphere_defaults['results_csv']
    )

    parser.add_argument(
        '--io_on',
        help="Enable or disable input/output. Accepts true/false, yes/no, 1/0.",
        type=parse_bool,
        default=dry_baroclinic_sphere_defaults['io_on']
    )

    args, unknown = parser.parse_known_args()

    dry_baroclinic_sphere(**vars(args))
