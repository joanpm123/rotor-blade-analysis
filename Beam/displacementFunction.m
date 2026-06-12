function ue = displacementFunction(data,Tn,u)

ue = zeros(data.ni*data.nne,data.nel);
I = zeros(data.ni*data.nne,1);

for e=1:data.nel
    for j=1:data.ni
        I(j,1) = nod2dof(data.ni,Tn(e,1),j);
        I(6+j,1) = nod2dof(data.ni,Tn(e,2),j);
    end
    ue(:,e) = u(I,1);
end
end