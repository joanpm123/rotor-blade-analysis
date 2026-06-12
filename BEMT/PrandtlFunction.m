function [F,inflow,phi,iterF] = PrandtlFunction(F,tolerance,H,sol,lmda_c,theta,rj,j)

er_F = 10;
iterF = 0;
max_iterF = 300;

while er_F>tolerance  && (iterF < max_iterF)
    iterF = iterF + 1;

    % Inflow calculation
    %{
    Aterm = H.sol*H.Cl_alpha/(16*F) - lmda_c/2;
    rad = Aterm^2 + H.sol*H.Cl_alpha/(8*F)*theta*rj(j);
    inflow = -Aterm + sqrt(rad);
    %}
    
    Aterm  = sol*H.Cl_alpha/(16*F);
    inflow = Aterm*(-1 + sqrt(1 + 2/Aterm*theta*rj(j)));

    %inflow = Aterm*(-1 + sqrt(1 + 32*F/(H.sol*H.Cl_alpha)*theta*rj(j)));

    % Induced inflow angle
    phi = inflow/rj(j);
    f   = H.Nb/2*(1-rj(j))/(rj(j)*phi);
    arg = exp(-f);

    % Prandtl's tip loss factor
    F_new = 2/pi*acos(arg);
    er_F = abs(F_new-F);
    F = F_new;
end

end