function [K,M] = assemblyStiffness(data,Td,Kel,Mel)

I=repmat(permute(Td,[2,3,1]),1,size(Td,2));
J=repmat(permute(Td,[3,2,1]),size(Td,2),1);

% Stiffness matrix
K=sparse(I(:),J(:),Kel(:),data.ndof,data.ndof);

% Mass matrix
M=sparse(I(:),J(:),Mel(:),data.ndof,data.ndof);
end