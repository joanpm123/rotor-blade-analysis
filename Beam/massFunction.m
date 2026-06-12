function [Mel,N] = massFunction(data,x,Tn,secProp,gp,gw,R)

Mel = zeros(data.nne*data.ni,data.nne*data.ni,data.nel);
N = zeros(data.ni,data.nne*data.ni,data.nel,size(gp,2));

% Unlike in shear stiffness calculation, here we use 2 Gauss points to
% provide exact integration

for e=1:data.nel

    l = norm(x(Tn(e,2),:)-x(Tn(e,1),:));     % Element size / Jacobian

    % Element properties
    rho = secProp.rho(e);
    A   = secProp.A(e);
    J   = secProp.J(e).*secProp.kt(e);
    Iy  = secProp.Iy(e);
    Iz  = secProp.Iz(e);
    
    r = rho*diag([A A A J Iy Iz]);

    % Identity matrix
    i = eye(6);
    
    % Loop over Gauss points
    for k=1:size(gp,2)                      
        % Shape functions at Gauss points
        N1 = (1-gp(1,k))/2;
        N2 = (1+gp(1,k))/2;
        N(:,:,e,k) = [N1*i N2*i];
        
        Mel(:,:,e) = Mel(:,:,e) + gw(1,k)*l*R(:,:,e)'*N(:,:,e,k)'*r*N(:,:,e,k)*R(:,:,e)/2;
    end
    
end

end