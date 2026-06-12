function [Ks,Bs] = shearFunction(data,l,R,dN,secProp,e)

% To avoid shear locking (induced by Timoshenko beam theory) we
% subintegrate the shear stiffness matrix, i.e. only using 1 Gauss
% integration point

% Shape function assuming only 1 Gauss point
N = 1/2;


% Strain matrix

Bs = sparse([1 1 1 1 2 2 2 2], ...
            [2 6 8 12 3 5 9 11], ...
            [dN(1) -N dN(2) -N dN(1) N dN(2) N],2,data.nne*data.ni);

% Element properties
G  = secProp.G(e);
A  = secProp.A(e);
ky = secProp.ky(e);
kz = secProp.kz(e);

% Constitutive matrix
Cs = G*A*[ky, 0; 0, kz];

% Shear component of stiffness matrix
Ks = l*R'*Bs'*Cs*Bs*R;
end