function blade = BladePlotAirfoil(displacement, propTable_airfoil, data, t,theta_x, H, theta_0, r_break, secProp, ea, eb, es, et)
% BLADEPLOTAIRFOIL  Render the deformed helicopter blade with airfoil
%                   variation and per-vertex stress, computed directly on
%                   the analytical NACA contour so geometry and stress
%                   are perfectly consistent.
%
%   displacement      : (data.nel x 3)
%   propTable_airfoil : 1x2 struct with field .naca_tc
%   data.nel          : number of beam elements
%   t                 : 1 -> sigma_xx,  else -> Von Mises
%   theta_x           : (data.nel x 1) elastic torsion [rad]
%   H                 : struct with H.theta_tw, H.c
%   theta_0           : collective pitch [rad]
%   r_break           : spanwise break station (non-dimensional)
%   secProp           : struct with .E(e), .G(e)
%   ea,eb,es,et       : strain arrays from strainFunction

    displacement = displacement(1:data.nel, :);
    theta_x      = theta_x(1:data.nel);
    nnod         = size(displacement, 1);

    % ---------------------------------------------------------------
    % 0. PITCH AND GEOMETRIC TWIST
    % ---------------------------------------------------------------
    theta_tw  = H.theta_tw * linspace(0, 1, nnod)';
    theta_def = -theta_0 - theta_tw + theta_x;

    % ---------------------------------------------------------------
    % 0b. CHORD PER BEAM NODE
    % ---------------------------------------------------------------
    c_vec = H.c(:);
    if numel(c_vec) ~= nnod
        c_vec = interp1(linspace(0,1,numel(c_vec)), c_vec, ...
                        linspace(0,1,nnod), 'linear')';
    end

    % ---------------------------------------------------------------
    % 1. THICKNESS RATIOS AND BLENDING
    % ---------------------------------------------------------------
    tc_root = propTable_airfoil(1).naca_tc;
    tc_tip  = propTable_airfoil(2).naca_tc;

    r_elem_vec = ((1:nnod)' - 0.5) / nnod;
    alpha = max(0, min(1, (r_elem_vec - r_break*0.8) ./ ...
                          (r_break*1.2 - r_break*0.8 + 1e-10)));

    % ---------------------------------------------------------------
    % 2. ANALYTICAL NACA CONTOUR
    %    N_prof points per surface, origin at c/4 (beam axis)
    % ---------------------------------------------------------------
    N_prof = 100;
    xi     = linspace(0, 1, N_prof)';

    a0 =  0.2969;  a1 = 0.1260;  a2 = 0.3516;
    a3 =  0.2843;  a4 = 0.1015;

    % Half-thickness at unit chord, scaled later per element
    t_unit = 5*(a0*sqrt(xi) - a1*xi - a2*xi.^2 + a3*xi.^3 - a4*xi.^4);

    nL = 2 * N_prof;   % contour points per cross-section

    % ---------------------------------------------------------------
    % 3. BUILD SKIN: geometry + stress evaluated AT contour points
    % ---------------------------------------------------------------
    skin_xyz   = zeros(nnod * nL, 3);
    skin_sigma = zeros(nnod * nL, 1);

    for i = 1:nnod

        % --- Blended section geometry ---
        a    = alpha(i);
        tc_i = (1-a)*tc_root + a*tc_tip;
        c_i  = c_vec(i);

        zt      = tc_i * c_i * t_unit;   % half-thickness [m]
        y_chord = xi * c_i - c_i/4;      % y ∈ [-c/4, 3c/4], origin at c/4

        % Closed contour: upper (LE→TE) then lower (TE→LE)
        y_contour = [y_chord;       flipud(y_chord)];   % [nL x 1]
        z_contour = [zt;           -flipud(zt)       ];

        % --- Stress evaluated at each contour point ---
        %     Same formulae as CrossSectionS2 — but now y,z are the
        %     exact contour coordinates, so stress and geometry match.
        E_e = secProp.E(i);
        G_e = secProp.G(i);

        y_e = y_contour;
        z_e = z_contour;

        exx = ea(1,i) + z_e*eb(1,i) - y_e*eb(2,i);
        gxy = es(1,i) - z_e*et(1,i);
        gxz = es(2,i) + y_e*et(1,i);

        sig_e = E_e * exx;
        txy_e = G_e * gxy;
        txz_e = G_e * gxz;
        vm_e  = sqrt(sig_e.^2 + 3*(txy_e.^2 + txz_e.^2));

        if t == 1
            stress_i = sig_e;
        else
            stress_i = vm_e;
        end

        % --- 3-D position (rotation by theta_def, translate to beam node) ---
        x0 = displacement(i, 1);
        y0 = displacement(i, 2);
        z0 = displacement(i, 3);
        ct = cos(theta_def(i));
        st = sin(theta_def(i));

        yr_d = ct*y_contour - st*z_contour;
        zr_d = st*y_contour + ct*z_contour;

        idx = (i-1)*nL + (1:nL);
        skin_xyz(idx, 1) = x0;
        skin_xyz(idx, 2) = y0 + yr_d;
        skin_xyz(idx, 3) = z0 + zr_d;
        skin_sigma(idx)  = stress_i;
    end

    % ---------------------------------------------------------------
    % 4. QUAD STRIP FACES
    % ---------------------------------------------------------------
    nFaces = (nnod-1) * nL;
    faces  = zeros(nFaces, 4);
    f = 1;
    for i = 1:(nnod-1)
        o1 = (i-1)*nL;
        o2 =  i   *nL;
        for j = 1:nL
            jn = mod(j, nL) + 1;
            faces(f,:) = [o1+j, o1+jn, o2+jn, o2+j];
            f = f + 1;
        end
    end

    % ---------------------------------------------------------------
    % 5. CAPS
    % ---------------------------------------------------------------
    root_face = 1:nL;
    tip_face  = (nnod-1)*nL + (1:nL);

    % ---------------------------------------------------------------
    % 6. PLOT
    % ---------------------------------------------------------------
    blade = figure('Color', 'w');

    patch('Faces',           faces, ...
          'Vertices',        skin_xyz, ...
          'FaceVertexCData', skin_sigma, ...
          'FaceColor',       'interp', ...
          'EdgeColor',       'none');
    hold on
    patch('Faces',           root_face, ...
          'Vertices',        skin_xyz, ...
          'FaceVertexCData', skin_sigma, ...
          'FaceColor',       'interp', ...
          'EdgeColor',       'none');
    patch('Faces',           tip_face, ...
          'Vertices',        skin_xyz, ...
          'FaceVertexCData', skin_sigma, ...
          'FaceColor',       'interp', ...
          'EdgeColor',       'none');
    hold off

    axis equal; view(3); grid on; box on

    xlabel('$X\,[m]$', 'Interpreter', 'latex', 'FontSize', 18)
    ylabel('$Y\,[m]$', 'Interpreter', 'latex', 'FontSize', 18)
    zlabel('$Z\,[m]$', 'Interpreter', 'latex', 'FontSize', 18)

    colormap(jet(256))
    cb = colorbar;
    cb.FontSize = 14;

    if t == 1
        ylabel(cb, '$\sigma_{xx}\,[Pa]$', 'Interpreter', 'latex', 'FontSize', 18)
        title('\textbf{Deformed Blade -- Airfoil Variation -- $\sigma_{xx}$ Distribution}', ...
              'Interpreter', 'latex', 'FontSize', 18)
    else
        ylabel(cb, '$\sigma_{VM}\,[Pa]$', 'Interpreter', 'latex', 'FontSize', 18)
        title('\textbf{Deformed Blade -- Airfoil Variation -- Von Mises Distribution}', ...
              'Interpreter', 'latex', 'FontSize', 18)
    end
end