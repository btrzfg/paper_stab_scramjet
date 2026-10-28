###################################################################################################
##                                                                                               ##
##        Hypersonic flow inside a 2-ramp scramjet model, two species, thermal equilibrium       ##
##                                                                                               ##
###################################################################################################

###################################################################################################
using OrdinaryDiffEqSSPRK
using Trixi
using AeroTrixi
using FlowRef
###################################################################################################
# Free-sream conditions:

const MASS_O2 = 5.3134e-26                  # O2 molecular mass [kg]
const Θ_VIBR_O2 = 2273.54                   # O2 characteristic vibrational temperature [K]
const X_O2 = 0.21                           # O2 molar fraction
const MASS_N2 = 4.65E-26                    # N2 molecular mass[ kg]
const Θ_VIBR_N2 = 3371                      # N2 characteristic vibrational temperature [K]
const X_N2 = 0.79                           # N2 molar fraction
const MASS_AVG = X_N2*MASS_N2+X_O2*MASS_O2  # mean molecular mass of the mixture [kg]
const Y_N2  = X_N2 * MASS_N2 / MASS_AVG     # N2 mass fraction [-]
const Y_O2  = X_O2 * MASS_O2 / MASS_AVG     # O2 mass fraction [-]
const V1_FREESTREAM = 2830.0                # m / s
const V2_FREESTREAM = 0.0
const P_FREESTREAM = 4100.0                 # Pa
const T_FREESTREAM = 370.0                  # K
const L_REF = 0.75                          # scramjet lenght [m]

const R_SPECIFIC = FlowRef.k_B / MASS_AVG
const RHO_FREESTREAM = P_FREESTREAM / (R_SPECIFIC * T_FREESTREAM)

# O.G.: every reference quantity follows from (p, T, rho, L); the equations are solved
#       in these units, so the free stream is O(1) by construction.
const REF_Q = p_T_rho_L(P_FREESTREAM, T_FREESTREAM, RHO_FREESTREAM, L_REF)

###################################################################################################
# Equations:
e_int_f_O2 = T -> e_rot_cont(MASS_O2, T) + e_vibr_iho(MASS_O2, Θ_VIBR_O2, T)
c_int_f_O2 = T -> c_rot_cont(MASS_O2, T) + c_vibr_iho(MASS_O2, Θ_VIBR_O2, T)

e_int_f_N2 = T -> e_rot_cont(MASS_N2, T) + e_vibr_iho(MASS_N2, Θ_VIBR_N2, T)
c_int_f_N2 = T -> c_rot_cont(MASS_N2, T) + c_vibr_iho(MASS_N2, Θ_VIBR_N2, T)

equations = CompressibleEulerEquationsMs1T2D(REF_Q, [MASS_O2, MASS_N2],
                                             [e_int_f_O2, e_int_f_N2], [c_int_f_O2, c_int_f_N2];
                                             T_min = 10.0, T_max = 4.0e4, dT = 1.0,
                                             min_T_jump = 1e-5)

###################################################################################################
# Initial conditions:

@inline function initial_condition_mach7_flow(x, t,
                                                   equations::CompressibleEulerEquationsMs1T2D)
    prim = SVector(V1_FREESTREAM / REF_Q.v_ref,
                   V2_FREESTREAM / REF_Q.v_ref,
                   T_FREESTREAM / REF_Q.T_ref,
                   Y_O2 * RHO_FREESTREAM / REF_Q.rho_ref,
                   Y_N2 * RHO_FREESTREAM / REF_Q.rho_ref)
    return prim2cons(prim, equations)
end

initial_condition = initial_condition_mach7_flow

@inline function boundary_condition_inflow(u_inner, normal_direction::AbstractVector, 
                                    x, t, surface_flux_function,
                                    equations::CompressibleEulerEquationsMs1T2D)
    
    u_boundary = initial_condition_mach7_flow(x, t, equations)
    flux = Trixi.flux(u_boundary, normal_direction, equations)
    return flux
end

@inline function boundary_condition_outflow(u_inner, normal_direction::AbstractVector, x, t, 
                                            surface_flux_function,
                                            equations::CompressibleEulerEquationsMs1T2D)
    
    flux = Trixi.flux(u_inner, normal_direction, equations)
    return flux
end

boundary_conditions = (
    PhysicalLine1 = boundary_condition_inflow, 
    PhysicalLine2 = boundary_condition_outflow, 
    PhysicalLine3 = boundary_condition_slip_wall, 
    PhysicalLine4 = boundary_condition_slip_wall)

###################################################################################################
#Spatial discretization:

polydeg = 2
surface_flux = FluxLaxFriedrichs(max_abs_speed)
volume_flux = flux_oblapenko
basis = LobattoLegendreBasis(polydeg)

shock_indicator = IndicatorHennemannGassner(equations, basis, alpha_max = 0.3, alpha_min = 0.001,
                                            alpha_smooth = true, variable = density_pressure)

volume_integral = VolumeIntegralShockCapturingHG(shock_indicator;
                                                 volume_flux_dg = volume_flux,
                                                 volume_flux_fv = surface_flux)

solver = DGSEM(polydeg = polydeg, surface_flux = surface_flux, volume_integral = volume_integral)

###################################################################################################
# MESH
# Generated in Gmsh and exported as .inp

mesh_file = joinpath(@__DIR__, "..", "..", "LorrainMesh_Trixi.inp")

mesh = P4estMesh{2}(mesh_file, polydeg = polydeg,
                    boundary_symbols=[:PhysicalLine1,:PhysicalLine2,:PhysicalLine3,:PhysicalLine4])

semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver,
                                    boundary_conditions = boundary_conditions)
###################################################################################################

function cons2prim_scaled(u, equations::CompressibleEulerEquationsMs1T2D)
    _, _, index_c, fracpos_c, prim = AeroTrixi.cons2prim_with_index(u, equations)
    v1, v2, T = prim[1], prim[2], prim[3]

    rho = density(u, equations)
    rho_inv = 1 / rho
    p = pressure(T, u, equations)

    gamma_T = AeroTrixi.get_gamma(u, rho_inv, index_c, fracpos_c, equations.thermodata)
    a = sqrt(gamma_T * p * rho_inv)

    return SVector(rho * REF_Q.rho_ref, v1 * REF_Q.v_ref, v2 * REF_Q.v_ref,
                   p * REF_Q.p_ref, T * REF_Q.T_ref,
                   sqrt(v1^2 + v2^2) / a, gamma_T,
                   u[4] * REF_Q.rho_ref, u[5] * REF_Q.rho_ref)
end

function Trixi.varnames(::typeof(cons2prim_scaled),
                        ::CompressibleEulerEquationsMs1T2D)
    return ("rho", "v1", "v2", "p", "T", "M", "gamma", "rho_O2", "rho_N2")
end

###################################################################################################                                                
# Callbacks

# (tmin,tmax)
tspan = (0.0, 0.2)

ode = semidiscretize(semi, tspan)

summary_callback = SummaryCallback()

analysis_callback = Trixi.AnalysisCallback(semi, interval = 100)

stepsize_callback = StepsizeCallback(cfl = 2.5)

output_directory = joinpath(@__DIR__, "out")

save_solution = SaveSolutionCallback(dt = 0.1,
                                     save_initial_solution = true,
                                     save_final_solution = true,
                                     solution_variables = cons2prim_scaled,
                                     output_directory = output_directory)

save_restart = SaveRestartCallback(interval = 0, save_final_restart = true,
                                   output_directory = joinpath(@__DIR__, "restart"))

callbacks = CallbackSet(summary_callback, analysis_callback, stepsize_callback, save_solution)

###################################################################################################
# Run the simulation:

sol = solve(ode, SSPRK53(); dt = 1.0,
            adaptive = false, 
            ode_default_options()...,
	    maxiters = 9999999,
            callback = callbacks);

###################################################################################################
#Print output:
summary_callback()
###################################################################################################
