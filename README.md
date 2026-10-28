# Hypersonic scramjet simulations with Trixi.jl

This repository contains the setup files to reproduce the 2-D scramjet test case of the paper
*"Efficient high-order Discontinuous Galerkin simulations of high-enthalpy flows"*
(B. Fagá, G. Oblapenko, Y. Zhu, D. Doehring, M. Torrilhon, RWTH Aachen).

The case is a Mach 7 inviscid air flow through a two-ramp scramjet inlet. It is computed with
two gas models:

- calorically perfect: single-species compressible Euler equations with constant ratio of specific heat,
  using the standard [Trixi.jl](https://github.com/trixi-framework/Trixi.jl) equations;
- thermally perfect: a two-species N2/O2 mixture in thermal equilibrium, with rotational and
  vibrational (infinite harmonic oscillator) internal energies. It uses the multi-species Euler
  equations `CompressibleEulerEquationsMs1T2D` and the entropy-conservative flux
  `flux_oblapenko` from [AeroTrixi.jl](https://github.com/ACoM-RWTH/AeroTrixi.jl).

## Usage

The setup was tested with Julia 1.12.

Clone the repository, start Julia from the root directory with `julia --project=.`, and run:

```julia
using Pkg
Pkg.instantiate()
```
Then run a case from each case directory:

```bash
cd cases/calorically_perfect/
julia --project=... scramjet_calorically_perfect.jl
```

and

```bash
cd cases/thermally_perfect/
julia --project=... cases/thermally_perfect/scramjet_thermally_perfect.jl
```

**NOTE**: Be aware that to run on a local machine, the cases may take several hours. So its recommended to use more threads if possible (check your system). In order to do that, simply add `--threads=N`. For example, if you are in root directory:

```bash
cd cases/calorically_perfect
julia --project=... --threads=4 scramjet_calorically_perfect.jl
```

The solution files are written to an `out` directory, and you can post-process them with
[Trixi2Vtk.jl](https://github.com/trixi-framework/Trixi2Vtk.jl) and view them in ParaView.

## Case description

### Free-stream conditions

| Quantity                | Value                         |
|-------------------------|-------------------------------|
| Freestream Velocity     | 2830 m/s                      |
| Freestream Pressure     | 4100 Pa                       |
| Freestream Temperature  | 370 K                         |
| Composition             | 79% N2, 21% O2 (molar)        |
| Reference length        | 0.75 m (scramjet length)      |

### Mesh

The unstructured quadrilateral mesh `LorrainMesh_Trixi.inp` was generated in Gmsh and its read
as a `P4estMesh{2}` with curved elements of the same degree as the solution.

### Numerical setup

Both cases use the same discretization:

- DGSEM with polynomial degree `polydeg = 2`;
- Lax–Friedrichs surface flux;
- shock capturing that blends DG and finite volumes, with the Hennemann–Gassner indicator on
  density pressure (`alpha_max = 0.3`, `alpha_min = 0.001`);
- SSPRK53 time integration with CFL = 2.5.

| | Calorically perfect | Thermally perfect |
|---|---|---|
| Equations | `CompressibleEulerEquations2D(1.4)` | `CompressibleEulerEquationsMs1T2D` (N2, O2) |
| Volume flux | `flux_ranocha` | `flux_oblapenko` |
| Units | dimensional (SI) | non-dimensional, scaled with free-stream conditions |
| Final time | 0.0008 s | 0.2 (non-dimensional) |

In the thermally perfect case, the internal energy and specific heat of each species are

```julia
e_int_f_O2 = T -> e_rot_cont(MASS_O2, T) + e_vibr_iho(MASS_O2, Θ_VIBR_O2, T)
c_int_f_O2 = T -> c_rot_cont(MASS_O2, T) + c_vibr_iho(MASS_O2, Θ_VIBR_O2, T)
```

with Θ_vibr = 2273.54 K for O2 and 3371 K for N2. They are tabulated for
10 K ≤ T ≤ 4 × 10⁴ K with ΔT = 1 K. The output contains the dimensional density, velocity,
pressure, temperature, Mach number, and the partial densities.

## Citing

If you use this code in your research, please cite:

```bibtex
@inproceedings{faga2026efficient,
  title     = {Efficient high-order {D}iscontinuous {G}alerkin simulations of high-enthalpy flows},
  author    = {Fagá, B. and Oblapenko, G. and Zhu, Y. and Doehring, D. and Torrilhon, M.},
  booktitle = {25th DGLR-Fachsymposium der STAB},
  address   = {Aachen, Germany},
  series    = {Notes on Numerical Fluid Mechanics and Multidisciplinary Design},
  publisher = {Springer},
  year      = {2026},
  note      = {To appear}
}
```
