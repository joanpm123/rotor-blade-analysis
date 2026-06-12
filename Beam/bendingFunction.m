function [Kb,Bb] = bendingFunction(data,l,R,dN,secProp,e)

% Strain matrix
Bb = sparse([1 1 2 2],[5 11 6 12],[dN(1) dN(2) dN(1) dN(2)],2,data.nne*data.ni);

% Element properties
E   = secProp.E(e);
Iy  = secProp.Iy(e);
Iz  = secProp.Iz(e);
Iyz = secProp.Iyz(e);

% Constitutive matrix
Cb = E*[Iy, Iyz; 
        Iyz, Iz];   % Exp.(3.11a)

% Bending component of stiffness matrix
Kb = l*R'*Bb'*Cb*Bb*R;

end