function theta_0 = initalBEMT(H,Den)

% Density
rho = Den.rho0*(1-Den.delta*Den.h/Den.T0).^(9.81/Den.R_air/Den.delta-1);

% Thrust requirements to Hover
T = H.W;                          % Thrust [N]
CT_req = T/(rho*H.A*H.Vtip^2);    % Required C_T

% Starting point
theta_0 = 6*CT_req/(H.sol(1)*H.Cl_alpha) - 3/4*H.theta_tw...
              + 3/2*sqrt(CT_req/2);
end