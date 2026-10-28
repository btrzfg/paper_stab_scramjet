###################################################################################################
##                                                                                               ##
##        Hypersonic flow inside a 2-ramp scramjet model, calorically perfect.                   ##
##                                                                                               ##
###################################################################################################

###################################################################################################
using Trixi
using OrdinaryDiffEqSSPRK
###################################################################################################
# Case conditions
equations = CompressibleEulerEquations2D(1.4)

const R_air = 287.0  # J/(kg·K), for the ideal gas law p = rho*R*T

# initial conditions for density, velocity, pressure
@inline function ic_mach7_flow(x, t, equations::CompressibleEulerEquations2D)
    rho_inf = 0.0385     # from p = rho*R*T ; T = 370
    v1      = 2830.0     # v_x m/s
    v2      = 0.0        # v_y m/s
    p_inf   = 4100.0     # Pa
    prim = SVector(rho_inf, v1, v2, p_inf)
    return prim2cons(prim, equations)
end

initial_condition = ic_mach7_flow

@inline function bc_inflow(u_inner, normal_direction::AbstractVector, x, t, surface_flux_function,
                                    equations::CompressibleEulerEquations2D)
    # inflow boundary condition... NO REFLECTION 
    u_boundary = ic_mach7_flow(x, t, equations)
    flux = Trixi.flux(u_boundary, normal_direction, equations)
    return flux
end

@inline function bc_outflow(u_inner, normal_direction::AbstractVector, x, t, surface_flux_function,
                                    equations::CompressibleEulerEquations2D)
    # outflow boundary condition... NO REFLECTION 
    flux = Trixi.flux(u_inner, normal_direction, equations)
    return flux
end

###################################################################################################
# Output variables written to disk.
# We return the usual primitives PLUS temperature, velocity magnitude and Mach.
@inline function cons2out(u, equations::CompressibleEulerEquations2D)
    rho, rho_v1, rho_v2, rho_e = u
    v1 = rho_v1 / rho
    v2 = rho_v2 / rho
    p  = Trixi.pressure(u, equations)

    v_mag = sqrt(v1^2 + v2^2)               # velocity magnitude [m/s]
    a     = sqrt(equations.gamma * p / rho) # local speed of sound [m/s]
    mach  = v_mag / a                       # local Mach number [-]
    T     = p / (rho * R_air)               # temperature [K], ideal gas law

    return SVector(rho, v1, v2, p, T, v_mag, mach)
end

Trixi.varnames(::typeof(cons2out), ::CompressibleEulerEquations2D) =
    ("rho", "v1", "v2", "p", "T", "v_mag", "Mach")

###################################################################################################
# spatial discretization
polydeg = 2
surface_flux = flux_lax_friedrichs 
volume_flux = flux_ranocha
basis = LobattoLegendreBasis(polydeg)
shock_indicator = IndicatorHennemannGassner(equations, basis, alpha_max = 0.3, alpha_min = 0.001, 
                                            alpha_smooth = true, variable = density_pressure)
# Look up the references to put here
volume_integral = VolumeIntegralShockCapturingHG(shock_indicator;
                                                 volume_flux_dg = volume_flux,
                                                 volume_flux_fv = surface_flux)
                                                 
###################################################################################################                                                
# Solver
solver = DGSEM(polydeg = polydeg, surface_flux = surface_flux, volume_integral = volume_integral)
###################################################################################################
#  mesh imported from Gmsh

mesh_file = joinpath(@__DIR__, "..", "..", "LorrainMesh_Trixi.inp")
 
mesh = P4estMesh{2}(mesh_file, polydeg = polydeg, 
                    boundary_symbols = [:PhysicalLine1, :PhysicalLine2, :PhysicalLine3, :PhysicalLine4])

bc_mach7_flow = (
    PhysicalLine1 = bc_inflow, 
    PhysicalLine2 = bc_outflow, 
    PhysicalLine3 = boundary_condition_slip_wall, 
    PhysicalLine4 = boundary_condition_slip_wall
)

semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver, 
                                    boundary_conditions = bc_mach7_flow)

tspan = (0.0, 0.0008 )
ode = semidiscretize(semi, tspan)

###################################################################################################
#callback functions during the time integration

summary_callback = SummaryCallback()

analysis_callback = AnalysisCallback(semi, interval = 100 )

stepsize_callback = StepsizeCallback(cfl = 2.5)

save_solution = SaveSolutionCallback(dt = 4.0E-4,
                                     solution_variables = cons2out) 
###################################################################################################
# Run the simulation

callbacks = CallbackSet(summary_callback, analysis_callback, stepsize_callback, save_solution)
sol = solve(ode, SSPRK53();
            dt = 1.0,
            ode_default_options()...,
            adaptive = false,
            callback = callbacks);

summary_callback()
###################################################################################################
