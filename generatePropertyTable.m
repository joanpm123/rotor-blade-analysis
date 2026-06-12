%% GENERATE CROSS-SECTION PROPERTY TABLE
% Calls the beam_section_properties mesher at multiple chord values
% (and optionally multiple NACA profiles) to build a lookup table.
%
% The beam axis lives at c/4, so all offsets (y_sc, y_cg) are stored
% relative to that line.
%
% OUTPUT: saves 'propertyTable.mat' containing struct array propTable
%         with fields per entry:
%   .chord    – chord length [m]
%   .naca_h   – NACA half-thickness parameter [m]
%   .A        – cross-section area [m^2]
%   .y_cg     – centroid y-offset from c/4 [m]
%   .z_cg     – centroid z-offset from c/4 [m]
%   .y_sc     – shear-centre y-offset from c/4 [m]
%   .z_sc     – shear-centre z-offset from c/4 [m]
%   .Iy_cm    – 2nd moment about y' through centroid [m^4]
%   .Iz_cm    – 2nd moment about z' through centroid [m^4]
%   .Iyz_cm   – product of inertia through centroid [m^4]
%   .Iy_c4    – 2nd moment about y' through c/4 [m^4]
%   .Iz_c4    – 2nd moment about z' through c/4 [m^4]
%   .Iyz_c4   – product of inertia through c/4 [m^4]
%   .J        – torsional constant = kt * Isc [m^4]
%   .Isc      – polar moment about SC [m^4]
%   .ky       – shear correction factor (chord direction)
%   .kz       – shear correction factor (thickness direction)
%   .kt       – torsion correction factor
%   .nodes    – cross-section mesh nodes (for stress recovery)
%   .elems    – cross-section mesh connectivity
%
% USAGE: Run this script from MATLAB with the mesher folder on the path.

clear; clc; close all;

% Add mesher to path
mesherPath = fullfile(fileparts(mfilename('fullpath')), ...
    '..', 'beam_section_properties_v2');
addpath(mesherPath);

%% ============ CONFIGURATION ============

% Conversion factor: mil to meters
miltom = 1.0*25.4e-5;

% --- Baseline section parameters (UH-60 Black Hawk) ---
c_base   = 0.527;        % Baseline chord [m]
h_ratio  = 0.145;        % NACA thickness-to-chord ratio (SC1095 ≈ 0.095, approx 0.145 here)
t_skin   = 40*miltom;    % Skin thickness [m] — FIXED (absolute)
t_spar1  = 200*miltom;   % Front spar thickness [m] — FIXED
t_spar2  = 200*miltom;   % Rear spar thickness [m] — FIXED
p1_ratio = 0.15;         % Front spar position (fraction of chord)
p2_ratio = 0.41;         % Rear spar position (fraction of chord)

% Material properties (same for all sections)
density      = 1551;
youngModulus = 50e9;
poissonRatio = 0.3;

% Mesh resolution
resolution = 1;

%% ============ STUDY 1: CHORD VARIATION ============
% Chord multipliers relative to baseline
chord_multipliers = [0.70, 0.85, 1.00, 1.15, 1.30];
chord_values = c_base * chord_multipliers;

fprintf('=== STUDY 1: Chord variation (fixed airfoil shape) ===\n\n');

nChord = length(chord_values);
propTable_chord = struct();

for i = 1:nChord
    c = chord_values(i);
    h = h_ratio * c;            % NACA half-thickness scales with chord
    p1 = p1_ratio * c;          % Spar positions scale with chord
    p2 = p2_ratio * c;
    % Skin and spar thicknesses remain FIXED (absolute)

    param = [c, h, t_skin, p1, t_spar1, p2, t_spar2];
    geometry = "spars";

    fprintf('  Computing chord = %.4f m (%.0f%% of baseline)...\n', ...
        c, chord_multipliers(i)*100);

    % Generate mesh
    [nodes, mat, elems, elemat] = getMeshData(geometry, param, ...
        resolution, density, youngModulus, poissonRatio, "", "", false);

    % Compute section properties
    section = getSectionProperties(nodes, mat, elems, elemat, false);

    % --- Store properties ---
    propTable_chord(i).chord  = c;
    propTable_chord(i).naca_h = h;

    % Area
    propTable_chord(i).A = section.A;

    % Centroid position from LE (mesher coordinates)
    y_cm_LE = section.xcm(1);    % chordwise, from LE
    z_cm_LE = section.xcm(2);

    % Shear centre position from LE
    y_sc_LE = section.xsc(1);
    z_sc_LE = section.xsc(2);

    % Offsets from sc (beam axis)
    y_c4 = c / 4;
    propTable_chord(i).y_c4 = y_c4 - y_sc_LE;
    propTable_chord(i).y_cg = y_cm_LE - y_sc_LE;
    propTable_chord(i).z_cg = z_cm_LE;           % = 0 for symmetric
    propTable_chord(i).y_sc = y_sc_LE;
    propTable_chord(i).z_sc = z_sc_LE;           % = 0 for symmetric

    % Inertias about centroid (mesher output)
    % Mesher convention: Ixx = ∫y² dA (about chord axis)  → Iy in beam
    %                    Iyy = ∫x² dA (about thickness axis) → Iz in beam
    propTable_chord(i).Iy_cm  = section.Ixx_0;
    propTable_chord(i).Iz_cm  = section.Iyy_0;
    propTable_chord(i).Iyz_cm = section.Ixy_0;

    % Inertias about sc via parallel axis theorem
    e_y = y_cm_LE - y_sc_LE;       % centroid offset from c/4 in chord direction
    e_z = z_cm_LE;                 % centroid offset from c/4 in thickness direction
    propTable_chord(i).Iy_sc  = section.Ixx_0 + e_z^2 * section.A;
    propTable_chord(i).Iz_sc  = section.Iyy_0 + e_y^2 * section.A;
    propTable_chord(i).Iyz_sc = section.Ixy_0 + e_y * e_z * section.A;

    % Torsion
    propTable_chord(i).J   = section.Isc;
    propTable_chord(i).kt  = section.kt;

    % Shear correction factors
    % Mesher kx → beam ky (chord dir), mesher ky → beam kz (thickness dir)
    propTable_chord(i).ky = section.kx;
    propTable_chord(i).kz = section.ky;

    % Cross-section mesh (for stress recovery)
    propTable_chord(i).nodes = nodes;
    propTable_chord(i).elems = elems;

    % Store raw LE-referenced positions for convenience
    propTable_chord(i).y_cm_LE = y_cm_LE;
    propTable_chord(i).y_sc_LE = y_sc_LE;

    fprintf('    A=%.4e  Iy=%.4e  Iz=%.4e  J=%.4e\n', ...
        section.A, section.Ixx_0, section.Iyy_0, section.kt*section.Isc);
    fprintf('    y_cg(sc)=%.4f  y_c4(sc)=%.4f  ky=%.4f  kz=%.4f  kt=%.5f\n', ...
        propTable_chord(i).y_cg, propTable_chord(i).y_c4, ...
        propTable_chord(i).ky, propTable_chord(i).kz, propTable_chord(i).kt);
end


%% ============ STUDY 2: AIRFOIL VARIATION ============
% Two NACA profiles: thick root, thin tip
% Candidates: NACA 0024 (root), NACA 0012 (tip)
% Both at baseline chord

fprintf('\n=== STUDY 2: Airfoil variation (fixed chord) ===\n\n');

naca_thicknesses = [0.18, 0.1];     % t/c ratios (NACA 00xx → xx/100)
naca_labels = {'NACA0024_root', 'NACA0012_tip'};

propTable_airfoil = struct();

for i = 1:length(naca_thicknesses)
    c = c_base;
    h = naca_thicknesses(i) * c;   % NACA half-thickness
    p1 = p1_ratio * c;
    p2 = p2_ratio * c;

    param = [c, h, t_skin, p1, t_spar1, p2, t_spar2];
    geometry = "spars";

    fprintf('  Computing %s (t/c = %.0f%%)...\n', naca_labels{i}, naca_thicknesses(i)*100);

    [nodes, mat, elems, elemat] = getMeshData(geometry, param, resolution, density, youngModulus, poissonRatio, "", "", false);
    section = getSectionProperties(nodes, mat, elems, elemat, false);

    y_cm_LE = section.xcm(1);
    z_cm_LE = section.xcm(2);
    y_sc_LE = section.xsc(1);
    z_sc_LE = section.xsc(2);
    y_c4 = c / 4;

    propTable_airfoil(i).chord    = c;
    propTable_airfoil(i).naca_tc  = naca_thicknesses(i);
    propTable_airfoil(i).naca_h   = h;
    propTable_airfoil(i).label    = naca_labels{i};
    propTable_airfoil(i).A        = section.A;
    propTable_airfoil(i).y_cg     = y_cm_LE - y_sc_LE;
    propTable_airfoil(i).z_cg     = z_cm_LE;
    propTable_airfoil(i).y_c4     = y_c4 - y_sc_LE;
    propTable_airfoil(i).z_c4     = z_sc_LE;
    propTable_airfoil(i).Iy_cm    = section.Ixx_0;
    propTable_airfoil(i).Iz_cm    = section.Iyy_0;
    propTable_airfoil(i).Iyz_cm   = section.Ixy_0;
    e_y = y_cm_LE - y_sc_LE;
    e_z = z_cm_LE;
    propTable_airfoil(i).Iy_sc    = section.Ixx_0 + e_z^2 * section.A;
    propTable_airfoil(i).Iz_sc    = section.Iyy_0 + e_y^2 * section.A;
    propTable_airfoil(i).Iyz_sc   = section.Ixy_0 + e_y*e_z * section.A;
    propTable_airfoil(i).J        = section.Isc;
    propTable_airfoil(i).kt       = section.kt;
    propTable_airfoil(i).ky       = section.kx;
    propTable_airfoil(i).kz       = section.ky;
    propTable_airfoil(i).nodes    = nodes;
    propTable_airfoil(i).elems    = elems;
    propTable_airfoil(i).y_cm_LE  = y_cm_LE;
    propTable_airfoil(i).y_sc_LE  = y_sc_LE;

    fprintf('    A=%.4e  Iy=%.4e  Iz=%.4e  J=%.4e\n', ...
        section.A, section.Ixx_0, section.Iyy_0, section.kt*section.Isc);
    fprintf('    y_cg(c/4)=%.4f  y_sc(c/4)=%.4f  ky=%.4f  kz=%.4f  kt=%.5f\n', ...
        propTable_airfoil(i).y_cg, propTable_airfoil(i).y_c4, ...
        propTable_airfoil(i).ky, propTable_airfoil(i).kz, propTable_airfoil(i).kt);
end


%% ============ SAVE ============
save('propertyTable.mat', 'propTable_chord', 'propTable_airfoil', ...
    'chord_values', 'chord_multipliers', 'naca_thicknesses', ...
    'c_base', 'h_ratio', 't_skin', 't_spar1', 't_spar2', ...
    'p1_ratio', 'p2_ratio', 'density', 'youngModulus', 'poissonRatio');

fprintf('\n  Property tables saved to propertyTable.mat\n');
