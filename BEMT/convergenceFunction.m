function [theta0_new,convergent] = convergenceFunction(H,error,tolerance,CT,CT_req,theta_0)

if error<=tolerance
    convergent = true;
    theta0_new = theta_0;
else
    convergent = false;

    theta0_new = theta_0 + 6*(CT_req-CT)/(H.sol(1)*H.Cl_alpha) +...
                3*sqrt(2)/4*(sqrt(CT_req) - sqrt(CT));
end
end