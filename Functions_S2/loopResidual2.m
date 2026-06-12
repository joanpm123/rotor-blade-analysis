function residual = loopResidual2(theta_y_vec, data, x, Tn, Td, Up, K, Mel, ...
        gp, gw, R, N_shape, pe, Be, H, Den, w, theta_x_in, theta_0,secProp,tolerance,study,y_cg_nod,y_c4_nod)
% LOOPRESIDUAL2  Residual for fsolve: converge BEMT-FEM at a given theta_y.

    theta_x     = theta_x_in;
    converg_loc = false;
    iter_loc    = 0;
    max_iter    = 50;
    theta_y     = theta_y_vec;   % initial guess 

    while ~converg_loc && (iter_loc < max_iter)

        theta_x_old = theta_x;

        [~, dT, converg_loc, theta0, ~, dD,~,~,~,~] = BEMT(data, H, Den, secProp, theta_0, theta_x, tolerance,study);
        theta_0 = theta0;
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

        Be_loc = Be;
        nBe = size(Be, 1);
        if iter_loc == 0
            Be_loc(nBe+1:nBe+2*data.nnod, :) = [(1:data.nnod)', 3*ones(data.nnod,1), w^2.*x(:,1).*theta_y;
                                                (1:data.nnod)', 4*ones(data.nnod,1), w^2.*x(:,1).*theta_y.*(y_cg_nod)];
        else
            Be_loc(nBe+1:nBe+2*data.nnod, :) = [(1:data.nnod)', 3*ones(data.nnod,1), w^2.*x(:,1).*theta_y_new;
                                                (1:data.nnod)', 4*ones(data.nnod,1), w^2.*x(:,1).*theta_y.*(y_cg_nod)];
        end

        % Element force vectors
        [Fe, Fel, ~, ~] = loadsFunction(data, x, Tn, pe, Qe, Be_loc, Mel, gw, gp, R, N_shape);

        % Loads, BC, solve
        f = assemblyLoads(data, Td, Fe, Fel);
        [ur, vr] = applyBC(data, Up);
        [u, ~]   = solveSystem(data, K, f, ur, vr);
        ue       = displacementFunction(data, Tn, u);
        
        % Update theta_x from torsion DOF
        theta_x = mean([ue(4,:); ue(10,:)], 1)';
        
        % Iteration index update
        iter_loc = iter_loc + 1;

        % Compute theta_y at all nodes from bending DOF
        theta_y_new = [ue(5,:)'; ue(11, data.nel)];

        % Secondary convergence: if theta_x stabilized, exit
        if iter_loc > 1 && norm(theta_x - theta_x_old) / (norm(theta_x_old) + 1e-15) < 1e-4
            break;
        end
    end
    residual = (theta_y_new(end) - theta_y_vec) / (abs(theta_y_vec) + 1e-15);
end