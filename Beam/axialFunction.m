function [Ka,Ba] = axialFunction(data,l,R,dN,secProp,e)

% Strain matrix
Ba = sparse(1,[1 7],[dN(1) dN(2)],1,data.nne*data.ni);

% Element material properties
E = secProp.E(e);
A = secProp.A(e);

% Constitutive matrix
Ca = E*A;

% Axial component of stiffness matrix
Ka = l*R'*Ba'*Ca*Ba*R;

end