function [u,r] = solveSystem(data,K,f,ur,vr)

% Solving the equation system using the partitioning method
vf    = setdiff(1:data.ndof,vr);              % Free DOFs array
u     = zeros(data.ndof,1);
u(vr) = ur;

% Free nodal displacements
u(vf) = K(vf,vf)\(f(vf)-K(vf,vr)*u(vr));      % Eq.(2.30) Belytschko

% Reaction loads at prescribed DOF
r     = K(vr,:)*u-f(vr);                      % Eq.(2.31) Belytschko
end