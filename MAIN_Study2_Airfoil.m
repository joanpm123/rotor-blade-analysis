%% PARAMETRIC STUDY 2: AIRFOIL VARIATION
% Unified BEMT-FEM analysis with beam axis at c/4.
% The blade has two spanwise zones with different NACA profiles,
% connected by a linear transition.
%
% Variables:
%   - Pair of NACA profiles (thick root / thin tip)
%   - Transition station r_break ∈ [0.4, 0.8]
%
% Held fixed: chord (baseline c_base), skin/spar thickness, material.
%
% PREREQUISITE: run generatePropertyTable.m to create propertyTable.mat

clear; clc; close all;

currentFolder = fileparts(mfilename('fullpath'));
addpath(genpath(fullfile(currentFolder, 'Beam')));
addpath(genpath(fullfile(currentFolder, 'BEMT')));
addpath(genpath(fullfile(currentFolder, 'Functions_S2')));

%% 1) PREPOCESS

%--------------------------- Input data ---------------------------------
load('propertyTable.mat');

% Helicopter and Blade parameters
HeliParameters

% Material Properties 
E   = youngModulus;                 % Young's modulus [Pa]
v   = poissonRatio;                 % Poisson's ratio
G   = E / (2*(1+v));                % Shear modulus [Pa]
rho = density;                      % Material density [kg/m^3]

%--------------------- Mesh Discretization ------------------------------
mesh_3D;

x = coordinates;

% Nodal connectivities
Tn = connectivity;

% Dimensions
data.nd   = size(x,2);              % Problem dimension
data.nel  = size(Tn,1);             % Number of elements
data.nnod = size(x,1);              % Number of nodes
data.nne  = size(Tn,2);             % Number of nodes in an element
data.ni   = size(x,2)*2;            % Degrees of freedom per node (3 dis-
                                    % placements and 3 rotations)
data.ndof = data.nnod*data.ni;      % Total number of degrees of freedom

% y' axis orientation
j = [0; 1; 0];


H.c   = c_base;
H.sol = (H.Nb.*H.c.*H.D/2)/H.A;     % Rotor Solidity
H.AR  = H.D/2./H.c;                 % Blade aspect ratio


%% ============ PARAMETRIC SWEEP CONFIGURATION ============

% Transition station (non-dimensional): where root airfoil ends and
% tip airfoil begins. Between r_break_start and r_break_end there is
% a linear transition.
r_break_values = [0.5 0.6 0.7 0.8 0.9];

% The two reference airfoil profiles come from propTable_airfoil:
%   propTable_airfoil(1) → thick root (e.g., NACA 0018)
%   propTable_airfoil(2) → thin tip   (e.g., NACA 0010)

study = 1;

%--------------------- Gauss quadrature definition ----------------------

gp = [-1/sqrt(3), 1/sqrt(3)];       % Gauss integration points
gw = [1,1];                         % Gauss weights


% Degrees of freedom connectivities matrix
Td = connectDOF(data,Tn);


%--------------------- Boundary Conditions ------------------------------

% Prescribed DOFs
Up = [1 1 0; 1 2 0; 1 3 0; 1 4 0; 1 5 0; 1 6 0];
                        % Global node, DOF, value
                        % Where, in the 2nd column, 1 refers to the dis-
                        % placement and 2 to the rotation. Displacement 
                        % DOFs are associated to point loads, while rota- 
                        % tions to bending moments.

data.np = size(Up,1);   % Number of prescribed nodes/DOFs


%% ============ PARAMETRIC LOOP ============

nCases  = length(r_break_values);
results = struct();


tolerance = 1e-6;


for iCase = 1:nCases

    fprintf('\n========== CASE %d / %d: r_break = %.2f ==========\n', ...
        iCase, nCases, r_break_values(iCase));

    r_break = r_break_values(iCase);

    % ======= Section Properties interpolation at element midpoint ======

    r_nodes = linspace(0, 1, nn)';
    r_elem  = (r_nodes(1:end-1) + r_nodes(2:end)) / 2;
    

    %---------------- Airfoil interpolation function --------------------
    P = interpAirfoil(propTable_airfoil, r_elem, r_break);

    secProp = struct();
    secProp.A      = P.A;
    secProp.E      = E * ones(data.nel, 1);
    secProp.G      = G * ones(data.nel, 1);
    secProp.rho    = rho * ones(data.nel, 1);
    secProp.Iy     = P.Iy_sc;
    secProp.Iz     = P.Iz_sc;
    secProp.Iyz    = P.Iyz_sc;
    secProp.J      = P.J;
    secProp.ky     = P.ky;
    secProp.kz     = P.kz;
    secProp.kt     = P.kt;
    secProp.y_c4   = P.y_c4;     % c/4 offset from SC [m]
    secProp.z_c4   = P.z_c4;
    secProp.y_cg   = P.y_cg;     % CG offset from SC [m]
    secProp.z_cg   = P.z_cg;
    secProp.chord  = H.c * ones(ne, 1);
    secProp.rbreak = r_break;

    % External point forces matrix: Global node, DOF, value [N]
    pe = [];                                      
    
    % Element stiffness matrix
    [Kel,Ka,Ba,Kb,Bb,Ks,Bs,Kt,Bt,R,T] = stiffnessFunction(data, x, Tn, secProp, j);

    % Element mass matrix
    [Mel, N] = massFunction(data, x, Tn, secProp, gp, gw, R);

    % Assemble global stiffness and mass matrices
    [K, M] = assemblyStiffness(data, Td, Kel, Mel);

    %% ----------- SOLVER (BEMT+FEM) -------------

    theta_x  = zeros(data.nel, 1);
    theta_0  = initalBEMT(H,Den);
    theta_y0 = -deg2rad(10);

    % Nodal CG offset
    y_cg_nod = zeros(data.nnod, 1);
    y_cg_nod(1:end-1) = y_cg_nod(1:end-1) + secProp.y_cg/2;
    y_cg_nod(2:end)   = y_cg_nod(2:end)   + secProp.y_cg/2;
    y_cg_nod(1,1)     = y_cg_nod(2,1); 

    % Nodal C/4 offset
    y_c4_nod = zeros(data.nnod, 1);
    y_c4_nod(1:end-1) = y_c4_nod(1:end-1) + secProp.y_c4/2;
    y_c4_nod(2:end)   = y_c4_nod(2:end)   + secProp.y_c4/2;
    y_c4_nod(1,1)     = y_c4_nod(2,1);
    

    % Body forces matrix: Global node, DOF, value [N/kg]
    w = H.Omega*2*pi/60;                % Rotor angular speed [rad/s]

    Be = [(1:data.nnod)', ones(data.nnod,1),  w^2 * x(:,1);
          (1:data.nnod)', 3*ones(data.nnod,1), -9.81 * ones(data.nnod,1);
          (1:data.nnod)', 4*ones(data.nnod,1), -9.81 * (y_cg_nod)];

    % --- fsolve for theta_y (loopResidual2 defined at end of file) ---
    options = optimoptions('fsolve', 'Display', 'iter', ...
        'TolFun', 1e-6, 'TolX', 1e-7, 'MaxIterations', 15);

    theta_y_sol = fsolve(@(ty) loopResidual2(ty, data, x, Tn, Td, Up, K, Mel, ...
        gp, gw, R, N, pe, Be, H, Den, w, theta_x, theta_0,secProp,tolerance,study,y_cg_nod,y_c4_nod), ...
        theta_y0, options);

    fprintf('Centrifugal projection angle: %.4f deg\n', rad2deg(theta_y_sol));

    %% --- Final solve ---

    theta_x = zeros(data.nel, 1);
    theta_0 = initalBEMT(H, Den);
    i_iter   = 0;
    iter_max = 50;
    converg  = false;
    nBe_base = size(Be, 1);

    while ~converg && (i_iter < iter_max)
        i_iter = i_iter + 1;

        Tni = zeros(data.nnod, 1);
        Dni = zeros(data.nnod, 1);

        [dP, dT, converg, theta0, CT, dD, dCP_i, dCP_0,~] = BEMT(data, H, Den, secProp, theta_0, theta_x, tolerance,study);
        theta_0 = theta0;
        dT = dT / H.Nb;

        Tni(1:end-1) = Tni(1:end-1) + dT/2;
        Tni(2:end)   = Tni(2:end)   + dT/2;
        Dni(1:end-1) = Dni(1:end-1) + dD/2;
        Dni(2:end)   = Dni(2:end)   + dD/2;

        % Distributed loads matrix    
        Qe = [(1:data.nnod)', 3*ones(data.nnod,1), Tni;
              (1:data.nnod)', 4*ones(data.nnod,1), Tni.*(-y_c4_nod);
              (1:data.nnod)', 2*ones(data.nnod,1), Dni];   % Global node, DOF, value [N/m]

        if i_iter == 1
            Be(nBe_base+1:nBe_base+2*data.nnod, :) = [(1:data.nnod)', 3*ones(data.nnod,1), w^2.*x(:,1).*theta_y_sol;
                                                      (1:data.nnod)', 4*ones(data.nnod,1), w^2.*x(:,1).*theta_y_sol.*(y_cg_nod)];
        else
            Be(nBe_base+1:nBe_base+2*data.nnod, :) = Mcf;
        end

        [Fe, Fel, ~, ~] = loadsFunction(data, x, Tn, pe, Qe, Be, Mel, gw, gp, R, N);
        f        = assemblyLoads(data, Td, Fe, Fel);
        [ur, vr] = applyBC(data, Up);
        [u, r_react] = solveSystem(data, K, f, ur, vr);
        ue       = displacementFunction(data, Tn, u);

        theta_x = mean([ue(4,:); ue(10,:)], 1)';
        theta_y = [ue(5,:)'; ue(11, data.nel)];

        %Mcf = [(1:data.nnod)', 3*ones(data.nnod,1), w^2*x(:,1).*theta_y];
        Mcf = [(1:data.nnod)', 3*ones(data.nnod,1), w^2.*x(:,1).*theta_y;
               (1:data.nnod)', 4*ones(data.nnod,1), w^2.*x(:,1).*theta_y_sol.*(y_cg_nod)];
    end

    % Figure of Merit
    CP_ideal = sum(dCP_i)/4;
    CP_0 = sum(dCP_0)/4;
    FM   = CP_ideal/(1.15*CP_ideal+CP_0);

    %% --- Store results ---

    results(iCase).r_break    = r_break;
    results(iCase).theta_0    = theta_0;
    results(iCase).theta_y    = theta_y(end);
    results(iCase).CT         = CT;
    results(iCase).u          = u;
    results(iCase).ue         = ue;
    results(iCase).theta_x    = theta_x;
    results(iCase).w_tip      = ue(3, end);
    results(iCase).twist_tip  = theta_x(end);
    results(iCase).secProp    = secProp;

    results(iCase).dCP_i      = dCP_i;
    results(iCase).dCP_0      = dCP_0;
    results(iCase).FM         = FM;

    [Fx, Fy, Fz, Mx, My, Mz] = forcesFunction(data,ue,Ka,Kb,Ks,Kt,R);
    results(iCase).Fx = Fx;
    results(iCase).Fy = Fy;
    results(iCase).Fz = Fz;
    results(iCase).Mx = Mx;
    results(iCase).My = My;
    results(iCase).Mz = Mz;

    [ea,eb,es,et] = strainFunction(data,ue,Ba,Bb,Bs,Bt,R);
    results(iCase).ea = ea;
    results(iCase).eb = eb;
    results(iCase).es = es;
    results(iCase).et = et;

    displacement = [ue(1:3,:)'; ue(7:9, data.nel)'] + x;
    results(iCase).displacement = displacement;

    fprintf('  theta_0 = %.4f deg, w_tip = %.4f m, twist_tip = %.4f deg\n', ...
        rad2deg(theta_0), ue(3,end), rad2deg(theta_x(end)));
end

%% ============ POSTPROCESS ============

fprintf('\n========== SUMMARY ==========\n');
fprintf('%12s %10s %10s %10s %12s\n', ...
    'r_break', 'theta_0', 'w_tip', 'twist_tip', 'CT');

for iCase = 1:nCases
    fprintf('%12.2f %10.4f %10.4f %10.4f %12.6f\n', ...
        results(iCase).r_break, rad2deg(results(iCase).theta_0), ...
        results(iCase).w_tip, rad2deg(results(iCase).twist_tip), ...
        results(iCase).CT);
end

%% --- Comparison plots ---

[Torsion_moment,Vertical_deflection,Power,FM,Angles,Axial] = ComparisonS2Function(results,nCases);

save('results_Study2.mat', 'results', 'data', 'x', 'Tn');
fprintf('\nResults saved to results_Study2.mat\n');


%% PLOTS FOR SPECIFIC BREAK POINTS

ratio = input(sprintf('Enter specific break point [%.2f, %.2f]: ', min(r_break_values), max(r_break_values)));

% Find row with that taper value
idx = find([results.r_break] == ratio);

% Check if it exists
if isempty(idx)
    error('Break point not found in results.')
end

ue           = results(idx).ue;
displacement = results(idx).displacement;
Fx           = results(idx).Fx;
Fy           = results(idx).Fy;
Fz           = results(idx).Fz;
Mx           = results(idx).Mx;
My           = results(idx).My;
Mz           = results(idx).Mz;
secProp      = results(idx).secProp;
ea           = results(idx).ea;
eb           = results(idx).eb;
es           = results(idx).es;
et           = results(idx).et;
theta_x      = results(idx).theta_x;
theta_0      = results(idx).theta_0;
H.c          = secProp.chord;


[sigma, tauy, tauz, VM] = CrossSectionS2(data, propTable_airfoil, secProp, ea, eb, es, et, r_break);

%% Cross-section stress distribution for a specific beam element (Study 2)

e = 1;   % <-- change to desired beam element

% --- Radial position of element e ---
r_elem_e = (e - 0.5) / data.nel;
r_break  = results(idx).r_break;

% --- Select mesh and connectivity for element e ---
if r_elem_e <= r_break
    nodes_e = propTable_airfoil(1).nodes;
    elems_e = propTable_airfoil(1).elems;
else
    nodes_e = propTable_airfoil(2).nodes;
    elems_e = propTable_airfoil(2).elems;
end

Ncs = size(nodes_e, 1);

% --- Stress fields (trim to actual mesh size) ---
sigma_xx = sigma(1:Ncs, 1, e);
tau_yz   = tauy(1:Ncs,  1, e);
tau_xz   = tauz(1:Ncs,  1, e);
VM_e     = VM(1:Ncs,    1, e);

% --- Plots ---
stress_data   = {sigma_xx,      tau_yz,       tau_xz,   VM_e};
stress_labels = {'\sigma_{xx}', '\tau_{xy}', '\tau_{xz}', '\sigma_{VM}'};

figure('Name', sprintf('Cross-section stresses — beam element %d', e))

for k = 1:4
    subplot(1, 4, k)
    patch('Faces', elems_e, 'Vertices', nodes_e, ...
          'FaceVertexCData', stress_data{k}, ...
          'FaceColor', 'interp', 'EdgeColor', 'none');
    axis equal
    colorbar
    xlabel('y'' [m]')
    ylabel('z'' [m]')
    title(sprintf('%s  —  element %d  (r = %.2f)', stress_labels{k}, e, r_elem_e))
end



%%
t=1;
blade = BladePlotAirfoil(displacement, propTable_airfoil, data, t,theta_x, H, theta_0, r_break, secProp, ea, eb, es, et);