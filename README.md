# BEMT + FEM Rotor Blade Analysis

**Joan Pratdepadua Mata** · Universitat Politècnica de Catalunya · June 2026

---

## 1. Overview

This repository contains a MATLAB implementation of a coupled **Blade Element Momentum Theory (BEMT)** and **Finite Element Method (FEM)** framework for the aeroelastic analysis of a helicopter rotor blade in hover. The structural model uses 150 two-node Timoshenko beam elements (6 DOF/node) with the beam axis at the shear centre. The aerodynamic model applies Prandtl tip-loss and iterates on the inflow ratio until convergence (tolerance 10⁻⁶). Aerodynamic loads (thrust, drag, pitching moment at c/4), centrifugal stiffening, and gravity are assembled into the FEM load vector at each aeroelastic iteration. The reference aircraft is the Sikorsky UH-60 Black Hawk (*D* = 16.36 m, 4 blades, 258 RPM).

Two parametric studies are provided. **Study 1** (`MAIN_Study1_Chord.m`) sweeps linear chord **taper ratios** from −0.20 to 0.35 (7 cases), examining how chord distribution affects rotor efficiency and blade structural response. **Study 2** (`MAIN_Study2_Airfoil.m`) varies the **spanwise airfoil transition** station *r*_break ∈ {0.5, 0.6, 0.7, 0.8, 0.9} between a thick root section (NACA 0018) and a thin tip section (NACA 0010).

---

## 2. Repository Structure

### Root files

| File | Description |
|------|-------------|
| `MAIN_Study1_Chord.m` | Entry point: chord-taper parametric study |
| `MAIN_Study2_Airfoil.m` | Entry point: airfoil-transition parametric study |
| `mesh_3D.m` | Timoshenko beam mesh generator (150 elements) |
| `generatePropertyTable.m` | Pre-processing: builds `propertyTable.mat` |
| `propertyTable.mat` | Pre-computed cross-section property tables |
| `nodes.mat` / `elems.mat` | 2-D cross-section FE mesh |
| `results_Study1.mat` / `results_Study2.mat` | Saved results |

### `BEMT/` — Aerodynamic module

| File | Description |
|------|-------------|
| `BEMT.m` | Core BEMT solver with Prandtl tip-loss |
| `HeliParameters.m` | Helicopter & air-density constants |
| `initalBEMT.m` | Initial collective pitch estimate |
| `convergenceFunction.m` | Inflow convergence loop |
| `PrandtlFunction.m` | Prandtl tip-loss factor |

### `Beam/` — Structural FEM module

| File | Description |
|------|-------------|
| `stiffnessFunction.m` | Element K matrix (axial, bending, shear, torsion) |
| `massFunction.m` | Consistent mass matrix (2-point Gauss) |
| `assemblyStiffness.m` | Global K & M assembly |
| `loadsFunction.m` | Distributed + body force integration |
| `solveSystem.m` / `applyBC.m` | BC application and linear solve |
| `strainFunction.m` / `forcesFunction.m` | Strain and internal force recovery |
| `figuresFunction.m` | FEM post-processing plots |

### `Functions_S1/` — Study 1 helpers (chord taper)

| File | Description |
|------|-------------|
| `interpolateProperties.m` | Chord-interpolated section properties |
| `loopResidual.m` | Residual function for θ_y fsolve |
| `CrossSection.m` | 2-D cross-section stress recovery |
| `BladePlotTaper.m` / `BladePlotDispTaper.m` | 3-D blade surface plots |

### `Functions_S2/` — Study 2 helpers (airfoil transition)

| File | Description |
|------|-------------|
| `interpAirfoil.m` | Spanwise airfoil interpolation |
| `loopResidual2.m` | Residual function for θ_y fsolve |
| `CrossSectionS2.m` | Stress recovery for mixed-airfoil blade |
| `BladePlotAirfoil.m` | 3-D blade surface plot |

---

## 3. How to Run

**Prerequisite:** MATLAB R2022b or later with the Optimization Toolbox (`fsolve`). All subfolders are added automatically via `addpath(genpath(...))`.

**Step 1 — Generate the property table (once).**
Run `generatePropertyTable.m` to compute cross-section properties (area, second moments, shear correction factors, torsion constant, CG and c/4 offsets) for all chord and airfoil variants and save them to `propertyTable.mat`.

```matlab
>> generatePropertyTable
```

**Step 2 — Study 1 (chord taper).**
Loops over 7 taper ratios, prints a summary table, and saves `results_Study1.mat`. After the loop the user is prompted to select a specific taper ratio for detailed plots.

```matlab
>> MAIN_Study1_Chord
```

**Step 3 — Study 2 (airfoil transition).**
Loops over 5 transition stations and saves `results_Study2.mat`. The user is then prompted for a specific *r*_break to display cross-section stress fields and 3-D blade plots.

```matlab
>> MAIN_Study2_Airfoil
```

---

## 4. Key Outputs

| Output | Description |
|--------|-------------|
| θ₀ | Collective pitch angle satisfying trimmed hover [rad] |
| w_tip | Flapwise tip deflection [m] |
| θ_x(r) | Spanwise elastic torsion distribution [rad] |
| C_T | Non-dimensional rotor thrust coefficient |
| FM | Figure of Merit = C_P,ideal / (1.15 C_P,ideal + C_P,0) |
| dT, dD | Spanwise thrust and drag per unit span [N/m] |
| Fx, Fy, Fz, Mx, My, Mz | Internal forces and moments along the beam [N], [N·m] |
| ε_a, ε_b, ε_s, ε_t | Axial, bending, shear, torsion strain components |
| σ_xx, τ_xy, τ_xz, σ_VM | Cross-section stress fields on 2-D FE mesh [Pa] |
| 3-D blade plots | Surface plots of stress or displacement over the deformed blade geometry |

