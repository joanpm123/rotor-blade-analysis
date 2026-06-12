function residual = loopResidual(theta_y_vec, data, x, Tn,Td, ~, Up, K, Mel, ...
                                  gp, gw, R, N, pe, Be, H, Den, w, secProp, theta_x_in, ...
                                  theta_0, iter,converg,tolerance,study,y_cg_nod,y_c4_nod)
theta_x = theta_x_in;
theta_y = theta_y_vec;   % current guess 

iter_max  = 15;

while converg == false && (iter < iter_max)


    % ----------- Blade Element Momentum Theory (BEMT) -------------
    [~,dT,converg,theta0,~,dD,~,~,~,~] = BEMT(data, H, Den, secProp, theta_0, theta_x, tolerance, study);
    
    % Root pitch angle update
    theta_0 = theta0;

    % Blade thrust gradient distribution
    dT = dT / H.Nb;

    % Element-wise nodal interpolation of Thrust and Drag
    Tni = zeros(data.nnod,1);       % Lift force nodal values matrix
    Dni = zeros(data.nnod,1);       % Drag force nodal values matrix

    Tni(1:end-1) = Tni(1:end-1) + dT/2;
    Tni(2:end)   = Tni(2:end)   + dT/2;

    Dni(1:end-1) = Dni(1:end-1) + dD/2;
    Dni(2:end)   = Dni(2:end)   + dD/2;

    % Distributed loads matrix
    Qe = [(1:data.nnod)', 3*ones(data.nnod,1), Tni;
          (1:data.nnod)', 4*ones(data.nnod,1), Tni.*(-y_c4_nod);
          (1:data.nnod)', 2*ones(data.nnod,1), Dni];

    if iter == 0
        filas = size(Be,1);
        Be(filas+1:filas+2*data.nnod, :) = [(1:data.nnod)', 3*ones(data.nnod,1),w^2.*x(:,1).*theta_y;
                                            (1:data.nnod)', 4*ones(data.nnod,1),w^2.*x(:,1).*theta_y.*(y_cg_nod)];
    else
        Be(filas+1:filas+2*data.nnod, :) = [(1:data.nnod)',3*ones(data.nnod,1),w^2*x(:,1).*theta_y_new;
                                            (1:data.nnod)',4*ones(data.nnod,1),w^2.*x(:,1).*theta_y.*(y_cg_nod)];
    end

    % Element force vectors
    [Fe, Fel, ~, ~] = loadsFunction(data, x, Tn, pe, Qe, Be, Mel, gw, gp, R, N);

    % Loads, BC, solve
    f        = assemblyLoads(data, Td, Fe, Fel);
    [ur, vr] = applyBC(data, Up);
    [u, ~]   = solveSystem(data, K, f, ur, vr);
    ue       = displacementFunction(data, Tn, u);

    % Update theta_x from torsion DOF
    theta_x = mean([ue(4,:); ue(10,:)], 1)';
    
    % Iteration index update
    iter = iter+1;

    % Compute theta_y at all nodes from bending DOF
    theta_y_new = [ue(5,:)' ; ue(11, data.nel)];   % (nnod x 1)
end
    % Residual: difference between new theta_y and current guess
    residual = (theta_y_new(end) - theta_y_vec)/abs(theta_y_vec  + 1e-15 );
end