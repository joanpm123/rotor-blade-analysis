function [Fx,Fy,Fz,Mx,My,Mz] = forcesFunction(data,ue,Ka,Kb,Ks,Kt,R)

Fx = zeros(1,data.nel);
Fy = zeros(1,data.nel);
Fz = zeros(1,data.nel);
Mx = zeros(1,data.nel);
My = zeros(1,data.nel);
Mz = zeros(1,data.nel);

for e=1:data.nel
    fint(:,e) = R(:,:,e)*(Ka(:,:,e)+Kb(:,:,e)+Ks(:,:,e)+Kt(:,:,e))*ue(:,e);
    Fx(:,e) = fint(1,e) - fint(7,e);
    Fy(:,e) = fint(2,e) - fint(8,e);
    Fz(:,e) = fint(3,e) - fint(9,e);
    Mx(:,e) = fint(4,e) - fint(10,e);
    My(:,e) = fint(5,e) - fint(11,e);
    Mz(:,e) = fint(6,e) - fint(12,e);
end
end