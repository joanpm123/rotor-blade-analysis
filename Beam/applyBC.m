function [ur,vr] = applyBC(data,Up)

ur=zeros(data.np,1);      % Displacement values vector
vr=zeros(data.np,1);      % DOF indexes of imposed displacements vector

for i=1:data.np
    vr(i)=nod2dof(data.ni,Up(i,1),Up(i,2));
    ur(i)=Up(i,3);
end
end