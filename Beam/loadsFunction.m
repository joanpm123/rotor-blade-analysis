function [Fe,Fel,Fvert,FvertT]=loadsFunction(data,x,Tn,pe,Qe,Be,Mel,gw,gp,R,N)

Fel = zeros(data.nne*data.ni,1,data.nel);

Fe = zeros(data.ndof,1);
Q  = zeros(data.nnod,data.ni);
B  = zeros(data.nnod,data.ni);

FvertT =zeros(1,data.nel);

% Point loads
for q=1:size(pe,1)
    I = nod2dof(data.ni,pe(q,1),pe(q,2));
    Fe(I,1) = Fe(I,1) + pe(q,3);
end

% Nodal distributed forces
for r=1:size(Qe,1)
    Q(Qe(r,1),Qe(r,2)) = Q(Qe(r,1),Qe(r,2))+Qe(r,3);
end

% Nodal body forces
for s=1:size(Be,1)
    B(Be(s,1),Be(s,2)) = B(Be(s,1),Be(s,2)) + Be(s,3);
end

for e=1:data.nel
    l = norm(x(Tn(e,2),:)-x(Tn(e,1),:));     % Element size

    b(:,e)  = [B(Tn(e,1),:),B(Tn(e,2),:)]';
    qe(:,e) = [Q(Tn(e,1),:),Q(Tn(e,2),:)]';
    
    Fel(:,1,e) = Mel(:,:,e)*b(:,e);
    
    % Vertical force associtated to body forces
    Fvert(:,e) = Fel(3,1,e)+Fel(9,1,e);                   

    % Loop over Gauss points
    for k=1:size(gp,2)
        a = gw(1,k)*l/2*R(:,:,e)'*N(:,:,e,k)'*N(:,:,e,k)*R(:,:,e)*qe(:,e);
        FvertT(:,e)=FvertT(:,e)+(a(3,1)+a(9,1));
        Fel(:,1,e) = Fel(:,1,e) + a;
    end

end
end