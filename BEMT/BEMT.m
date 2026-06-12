function [dCP,dT,convergent,theta0_new,fig2,dD,dCP_i,dCP_0,CT,theta_r] = BEMT(data,H,Den,secProp,theta_0,theta_x,tolerance,study)
%% DATA SETS & PROBLEM DEFINITON

% Density
rho = Den.rho0*(1-Den.delta*Den.h/Den.T0).^(9.81/Den.R_air/Den.delta-1);

% Thrust requirements to Hover
T = H.W;                               % Thrust [N]
CT_req = T/(rho*H.A*H.Vtip^2);         % Required C_T

% Blade discretization
ne      = data.nel;                    % Number of elements
r       = linspace(0,1,ne+1);
rj      = zeros(ne,1);                 % Centre position of each element
delta_r = zeros(ne,1);

% Values for the climb inflow ratio [vc = 0 bc. HOVER FLIGHT]   
vc = 0;
lmda_c = vc/H.Vtip; 

% Result vectors initialization
dCT   = zeros(ne,1);
dCP   = zeros(ne,1);
dCP_i = zeros(ne,1);
dCP_0 = zeros(ne,1);
dD    = zeros(ne,1);

max_iterF = 300;


%% COMPUTATION

lmda_ih = sqrt(CT_req/2);
lc = lmda_c/lmda_ih;              % Adimensional axial inflow ratio
    
% Adim. required induced inflow ratio 
if lc>=0
    li = -lc/2 + sqrt((lc^2)/4+1); 
elseif lc<=-2
    li = -lc/2 - sqrt((lc^2)/4-1);
end

CT = 0; CP = 0;

% Loop over elements
for j = 1:ne 

    delta_r(j) = r(j+1) - r(j);

    % Element centre position
    rj(j) = delta_r(j)/2 + r(j);
    
    % Element chord / solidity / aspect ratio / parasite drag coefficient
    if study == 1
        c   = H.c;
        sol = H.sol;
        AR  = H.AR;
        if rj(j)<=secProp.rbreak
            Cd0 = H.naca18;
        else
            Cd0 = H.naca10;
        end
    else
        c   = secProp.chord(j);
        sol = mean([H.sol(j) H.sol(j+1)]);
        AR  = mean([H.AR(j) H.AR(j+1)]);
        Cd0 = H.Cd0;
    end

    % Element pitch angle
    theta = theta_0 + rj(j)*H.theta_tw + theta_x(j);

    % Initializing Prandtl's correction factor
    F = 1;

    % Prandtl's tip loss correction factor implementation
    [F,inflow,phi,iterF] = PrandtlFunction(F,tolerance,H,sol,lmda_c,theta,rj,j);

    % Induced inflow ratio
    lmda_i = inflow - lmda_c;
    
    % If !convergence on F
    if iterF >= max_iterF
        warning('Element %d: Prandtl iteration did not converge (case %d). Using last F value.', j, 1);
    end
    
    % Element thrust and power coefficients
    dCT(j,1)   = 4*F*lmda_i*inflow*rj(j);
    dCP_i(j,1) = rho*H.A*H.Vtip^3*lmda_i*dCT(j,1)*delta_r(j)/H.Nb;
    dCP_c      = lmda_c*dCT(j,1);
    dCP_0(j,1) = rho*H.A*H.Vtip^3*sol/2*Cd0*r(j)^3*delta_r(j)/H.Nb; %*secProp.A(j)/c;
    dCP(j,1)   = dCP_i(j,1) + dCP_c + dCP_0(j,1);
    
    % Drag force [Exp.66 Helicopters Notes]
    Cl      = H.Cl_alpha*(theta - phi - H.alph0);
    k       = 1/(pi*AR*H.e);
    Cd      = H.Cd0 + k*Cl^2;
    dD(j,1) = 1/2*rho*(H.Omega*2*pi/60*rj(j)*(H.D/2))^2*c*Cd; % dD/dx

    % Total blade thrust and power coefficients
    CT = CT + dCT(j,1)*delta_r(j);
    CP = CP + dCP(j,1)*delta_r(j);

    theta_r(j,1) = theta;
end


dT1    = rho*H.A*H.Vtip^2*dCT;        % This is equivalent to dT/dr

dT     = dT1/(H.D/2);                 % This is equivalent to dT/dx


%% CONVERGENCE VERIFICATION 

error = abs(CT - CT_req);

[theta0_new,convergent] = convergenceFunction(H,error,tolerance,CT,CT_req,theta_0);


%% PLOTS
fig2 = [];
end
