function [sigma,tauy,tauz,VM] = CrossSection(data,nodes,secProp,ea,eb,es,et)

sigma = zeros(size(nodes,1),1,data.nel);
tauy  = zeros(size(nodes,1),1,data.nel);
tauz  = zeros(size(nodes,1),1,data.nel);
VM = zeros(size(nodes,1),1,data.nel);

for e = 1:data.nel
    
    % Cross section deformations
    exx = ea(1,e) + nodes(:,2)*eb(1,e) - nodes(:,1)*eb(2,e);
    gxy = es(1,e) - nodes(:,2)*et(1,e);
    gxz = es(2,e) + nodes(:,1)*et(1,e);

    % Cross section material properties
    E = secProp.E(e);                           % Young's modulus
    G = secProp.G(e);                           % Shear modulus

    % Cross section stresses
    sigma(:,1,e) = E*exx;
    tauy(:,1,e)  = G*gxy;
    tauz(:,1,e)  = G*gxz;

    % Von Misses Criterion
    VM(:,1,e) = sqrt(sigma(:,1,e).^2 + 3*(tauy(:,1,e).^2 + tauz(:,1,e).^2));
end

end

