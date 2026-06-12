%% ANALYSIS OF A ROTOR BLADE USING FEM AND BEMT 
% Algorithm for solving 3D Aero-beam coupling problem for a helicopter
% rotor blade
%
% The beam axis is located at the shear centre sc
%
% All bibliographical references, unless specified, are retrieved from 
% Structural Analysis with the FEM, Eugenio Oñate, Vol.2

clear,clc
close all

currentFolder = fileparts(mfilename('fullpath'));
addpath(genpath(fullfile(currentFolder, 'Beam')));
addpath(genpath(fullfile(currentFolder, 'BEMT')));
addpath(genpath(fullfile(currentFolder, 'Functions_S1')));

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

% ======= TAPER: linear taper ratio ∈ [0, 0.3] =======
taper = [-0.2 -0.1 0 0.1 0.2 0.3 0.35];

r0    = 0;
r_nod = linspace(0,1,data.nnod);

%% ============ PARAMETRIC LOOP ============

nCases  = length(taper);
results = struct();

for iCase = 1:nCases

    fprintf('\n========== CASE %d / %d: r_break = %.2f ==========\n',iCase, nCases, taper(iCase));
    
    tap = taper(iCase);

    % Chord distribution
    H.c   = 1*c_base*(1 - tap * (r_nod - r0)/(1 - r0)); 
    H.sol = (H.Nb.*H.c.*H.D/2)/H.A;                         % Rotor Solidity
    H.AR  = H.D/2./H.c;                                     % Blade aspect ratio

    fprintf('Linear taper: taper = %.2f, c_root = %.4f, c_tip = %.4f\n',tap, H.c(1), H.c(end));

    % ======= Section Properties interpolation at element midpoint ======

    % Element chord (average of node chords)
    c_el = (H.c(1:end-1) + H.c(2:end)) / 2;

    P = interpolateProperties(propTable_chord, c_el, 'chord');

    % Struct for element section properties
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
    secProp.z_c4   = P.z_sc;
    secProp.y_cg   = P.y_cg;     % CG offset from SC [m]
    secProp.z_cg   = P.z_cg;
    secProp.chord  = c_el;

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

%--------------------- External Forces ----------------------------------

    % External point forces matrix: Global node, DOF, value [N]
    pe = [];                                      

    % Body forces matrix: Global node, DOF, value [N/kg]
    w = H.Omega*2*pi/60;                % Rotor angular speed [rad/s]

    % Nodal CG offset
    y_cg_nod = zeros(data.nnod, 1);
    y_cg_nod(1:end-1) = y_cg_nod(1:end-1) + secProp.y_cg'/2;
    y_cg_nod(2:end)   = y_cg_nod(2:end)   + secProp.y_cg'/2;
    y_cg_nod(1,1)     = y_cg_nod(2,1); 

    % Nodal C/4 offset
    y_c4_nod = zeros(data.nnod, 1);
    y_c4_nod(1:end-1) = y_c4_nod(1:end-1) + secProp.y_c4'/2;
    y_c4_nod(2:end)   = y_c4_nod(2:end)   + secProp.y_c4'/2;
    y_c4_nod(1,1)     = y_c4_nod(2,1); 


    Be = [(1:data.nnod)',ones(data.nnod,1),w^2*x(:,1);   % Centrifugal force Omega^2*R 
          (1:data.nnod)',3*ones(data.nnod,1),-9.81*ones(data.nnod,1);
          (1:data.nnod)',4*ones(data.nnod,1),-9.81*ones(data.nnod,1).*(y_cg_nod)]; % Gravity 


    %% 2) SOLVER
    tic
    disp(size(Tn))

    % Element stiffness matrix
    [Kel,Ka,Ba,Kb,Bb,Ks,Bs,Kt,Bt,R,T] = stiffnessFunction(data,x,Tn,secProp,j);

    % Element mass matrix
    [Mel,N] = massFunction(data,x,Tn,secProp,gp,gw,R);

    % Assemble global stiffness and mass matrices
    [K,M] = assemblyStiffness(data,Td,Kel,Mel);

    % ----------- Blade Element Momentum Theory (BEMT) -------------

    theta_x   = zeros(data.nel,1);     % Initial torsion pitch angle
    theta_0   = initalBEMT(H,Den);     % Initial theta0 value

    % Initial guess for theta_y (for fsolve)
    theta_y0  = -deg2rad(4);          

    % Convergence parameters
    converg   = false;
    tolerance = 1e-6;
    iter      = 0;   
    study     = 2;

    % Helicopter Blade mass [kg]
    blm = 90;

    % ----------------- fsolve options --------------------
    options = optimoptions('fsolve','Display','iter','TolFun',1e-7,'TolX',1e-7, ...
                       'MaxIterations',  30);

    theta_y_sol = fsolve( ...
        @(theta_y_vec) loopResidual(theta_y_vec, data, x, Tn,Td, blm, Up, K, Mel,gp, gw, R, N, pe, Be, H,Den, w, ...
                                secProp, theta_x, theta_0, iter,converg,tolerance,study,y_cg_nod,y_c4_nod),theta_y0, options);


    fprintf('The centrifugal force vertical contribution projection angle %.4f deg.\n',rad2deg(theta_y_sol));


%% 2.2) SOLVER FOR PARAMETERS USING theta_y_sol AS PROJECTION ANGLE FOR Fcf
% ----------- Blade Element Momentum Theory (BEMT) -------------

    theta_x = zeros(data.nel,1);        % Initial torsion pitch angle
    theta_y = zeros(data.nel,1);        % Initial bending angle
    theta_0 = initalBEMT(H,Den);        % Initial theta0 value

    % Iteration parameters
    i = 0;
    iter_max  = 15;
    tolerance = 10^-6;                  % Analysis tolerance
    converg   = false;

    while converg == false && (i < iter_max)
        i = i+1;                        % Iteration index

        Tni = zeros(data.nnod,1);       % Lift force nodal values matrix
        Dni = zeros(data.nnod,1);       % Drag force nodal values matrix

        % Blade element momentum theory computation
        [dP,dT,converg,theta0,CT_dis,dD,dCP_i,dCP_0,CT,theta_r] = BEMT(data,H,Den,secProp,theta_0,theta_x,tolerance,study);
    
        % Root pitch angle update
        theta_0 = theta0;

        % Blade thrust gradient distribution
        dT = dT/H.Nb;
% ---------------------------------------------------------------
        % Element wise interpolation of thrust and drag
        Tni(1:end-1) = Tni(1:end-1) + dT/2;
        Tni(2:end)   = Tni(2:end)   + dT/2;

        Dni(1:end-1) = Dni(1:end-1) + dD/2;
        Dni(2:end)   = Dni(2:end)   + dD/2;

        % Distributed loads matrix    
        Qe = [(1:data.nnod)', 3*ones(data.nnod,1), Tni;
              (1:data.nnod)', 4*ones(data.nnod,1), Tni.*(-y_c4_nod);
              (1:data.nnod)', 2*ones(data.nnod,1), Dni];   % Global node, DOF, value [N/m]
    
        if i==1
            filas = size(Be,1);
            Be(filas+1:filas+2*data.nnod, :) = [(1:data.nnod)',3*ones(data.nnod,1),w^2.*x(:,1).*theta_y_sol;
                                                (1:data.nnod)',4*ones(data.nnod,1),w^2.*x(:,1).*theta_y_sol.*(y_cg_nod)];
        else
            Be(filas+1:filas+2*data.nnod, :) = Mcf;
        end
    
        % Prescribed loads
        [Fe,Fel,Fvert,FvertT] = loadsFunction(data,x,Tn,pe,Qe,Be,Mel,gw,gp,R,N);

        % Assemble loads
        f = assemblyLoads(data,Td,Fe,Fel);

        % Apply prescribed DOFs
        [ur,vr] = applyBC(data,Up);

        % Solve system
        [u,r] = solveSystem(data,K,f,ur,vr);

        % Element's displacements and rotations
        ue = displacementFunction(data,Tn,u);
    
        % Pitch angle contribution due to blade torsion    
        theta_x = mean([ue(4,:); ue(10,:)], 1)'; % Rotation values about the neutral axis (x axis)
    
        % Bending angle
        theta_y = [ue(5,:)' ; ue(11,data.nel)];
    
        % Blade vertical deflection
        def = [ue(3,:)' ; ue(9,data.nel)];

        tip(i)   = ue(3,end);
        angle(i) = ue(5,end);

        % Centrifugal stiffening (vertical contribution due to bending)
        Mcf = [(1:data.nnod)',3*ones(data.nnod,1),w^2.*x(:,1).*theta_y;
               (1:data.nnod)',4*ones(data.nnod,1),w^2.*x(:,1).*theta_y_sol.*(y_cg_nod)]; %*blm/(H.D/2)
    
    end

    toc

    %% 3) POSTPROCESS

    % Strain component
    [ea,eb,es,et] = strainFunction(data,ue,Ba,Bb,Bs,Bt,R);

    % Internal forces and moments
    [Fx,Fy,Fz,Mx,My,Mz] = forcesFunction(data,ue,Ka,Kb,Ks,Kt,R);

    displacement = [ue(1:3,:)' ; ue(7:9,data.nel)'] + x;

    
    % Figure of Merit
    CP_ideal = sum(dCP_i)/4;
    CP_0 = sum(dCP_0)/4;
    FM   = CP_ideal/(1.15*CP_ideal+CP_0);

    %% --- Store results ---

    results(iCase).taper      = tap;
    results(iCase).theta_0    = theta_0;
    results(iCase).theta_y    = theta_y(end);
    results(iCase).dT         = dT;
    results(iCase).dD         = dD;
    results(iCase).CT         = CT;
    results(iCase).theta_r    = theta_r;
    results(iCase).u          = u;
    results(iCase).ue         = ue;
    results(iCase).theta_x    = theta_x;
    results(iCase).w_tip      = ue(3, end);
    results(iCase).twist_tip  = theta_x(end);
    results(iCase).secProp    = secProp;
    
    results(iCase).dCP_i      = dCP_i;
    results(iCase).dCP_0      = dCP_0;
    results(iCase).FM         = FM;

    results(iCase).Fx = Fx;
    results(iCase).Fy = Fy;
    results(iCase).Fz = Fz;
    results(iCase).Mx = Mx;
    results(iCase).My = My;
    results(iCase).Mz = Mz;

    results(iCase).ea = ea;
    results(iCase).eb = eb;
    results(iCase).es = es;
    results(iCase).et = et;

    results(iCase).displacement = displacement;

    fprintf('  theta_0 = %.4f deg, w_tip = %.4f m, twist_tip = %.4f deg\n', ...
        rad2deg(theta_0), ue(3,end), rad2deg(theta_x(end)));
end

fprintf('\n========== SUMMARY ==========\n');
fprintf('%12s %10s %10s %10s %10s %10s\n','taper', 'theta_0', 'w_tip', 'twist_tip', 'CT', 'FM');

for iCase = 1:nCases
    fprintf('%12.2f %10.4f %10.4f %10.4f %10.4f %10.4f\n',results(iCase).taper, rad2deg(results(iCase).theta_0), ...
        results(iCase).w_tip, rad2deg(results(iCase).twist_tip), results(iCase).CT, results(iCase).FM);
end

%% --- Comparison plots ---

[Torsion_moment,Vertical_deflection,Power,FM,Angles,Axial] = ComparisonS1Function(results,nCases);

save('results_Study1.mat', 'results', 'data', 'x', 'Tn');
fprintf('\nResults saved to results_Study1.mat\n');

%% PLOTS FOR SPECIFIC TAPER RATIOS

ratio = input(sprintf('Enter specific taper ratio [%.2f, %.2f]: ', min(taper), max(taper)));

% Find row with that taper value
idx = find([results.taper] == ratio);

% Check if it exists
if isempty(idx)
    error('Taper ratio not found in results.')
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

% Internal forces and moments distribution
[figure1,figure2,figure3,figure4] = figuresFunction(x,ue,displacement,Fx,Fy,Fz,Mx,My,Mz);


%% AERODYNAMIC RESULTS
figure
plot(x(1:end-1,1),results(idx).dT*H.D/(2*data.nnod),'LineWidth',2)
xlim([x(1,1), x(end-1,1)]);
set(gca, 'FontSize', 12);
title('$\textbf{Lift spanwise distribution}$', 'Interpreter', 'latex', 'FontSize', 18)
xlabel('Span x [m]', 'Interpreter', 'latex', 'FontSize', 16.5)
ylabel('dL [N]', 'Interpreter', 'latex', 'FontSize', 16.5)
grid on
figure
plot(x(1:end-1,1),results(idx).dD*H.D/(2*data.nnod),'LineWidth',2)
xlim([x(1,1), x(end-1,1)]);
set(gca, 'FontSize', 12);
title('$\textbf{Drag spanwise distribution}$', 'Interpreter', 'latex', 'FontSize', 18)
xlabel('Span x [m]', 'Interpreter', 'latex', 'FontSize', 16.5)
ylabel('dD [N]', 'Interpreter', 'latex', 'FontSize', 16.5)
grid on

%% 3.1) CROSS-SECTION STRAIN AND STRESS

load('elems.mat');
load('nodes.mat');

% Elementc cross-section normal stress
[sigma,tauy,tauz,VM] = CrossSection(data,nodes,secProp,ea,eb,es,et);

%% Element 'e' Cross section sigma_xx / tau_yz / VM distribution
e        = 1; 
radial   = x(e,1);
sigma_xx = sigma(:,1,e);
tau_xy   = tauy(:,1,e);
tau_xz   = tauz(:,1,e);
VM_e     = VM(:,1,e);

figure
patch('Faces', elems,'Vertices', nodes,'FaceVertexCData', sigma_xx, ...
      'FaceColor','interp','EdgeColor','none');
colormap(jet(256))   

ylim([-0.1 0.1])
daspect([1 1 1])
set(gca, 'FontSize', 10);
cb = colorbar;      

ylabel(cb, 'Stress [Pa]', 'Interpreter', 'latex', 'FontSize', 11)
title(sprintf('$\\mathbf{\\sigma_{xx}\\ distribution\\ at\\ x = %.1f}$', radial), ...
    'Interpreter','latex','FontSize',11)
xlabel('y', 'Interpreter', 'latex', 'FontSize', 11)
ylabel('z', 'Interpreter', 'latex', 'FontSize', 11)

grid on

%% 3D Blade

k = 1;

if k==1
    % Sigma_xx distribution
    blade = BladePlotTaper(displacement, sigma, elems, nodes, data, k, theta_x, H, theta_0);
else
    % Von Mises distribution
    blade = BladePlotTaper(displacement, VM, elems, nodes, data, k, theta_x, H, theta_0);
end

% Displacement distribution
blade1 = BladePlotDispTaper(displacement, theta_x, elems, nodes, data, H, theta_0);
