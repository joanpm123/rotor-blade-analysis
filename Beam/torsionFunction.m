function [Kt,Bt] = torsionFunction(data,l,R,dN,secProp,e)

% Strain matrix
Bt = sparse(1,[4 10],[dN(1) dN(2)],1,data.nne*data.ni);


% Element properties
G  = secProp.G(e);
J  = secProp.J(e);
kt = secProp.kt(e);

% Constitutive value
Ct = G*J*kt;

% Torsion component of stiffness matrix
Kt = l*R'*Bt'*Ct*Bt*R;

end