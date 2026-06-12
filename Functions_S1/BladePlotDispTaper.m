function blade = BladePlotDispTaper(displacement, theta_x, elems, nodes, data, H, theta_0)
% BLADEPLOTDISPTAPER  Render the deformed TAPERED helicopter blade colored
%                     by the vertical (Z) displacement of every skin vertex.
%
%   displacement : (data.nel x 3)   [x y z] final positions of each beam node
%   theta_x      : (data.nel x 1)   elastic torsion angle [rad] at each beam node
%   elems        : (ne x k)         cross-section connectivity (zero-padded)
%                                   of the ROOT cross-section only
%   nodes        : (ncs x 2)        ROOT cross-section [y z] coordinates
%   data.nel     : number of beam nodes
%   H            : struct containing
%                    H.theta_tw  -> tip geometric twist [rad] (scalar; the
%                                   twist is distributed linearly from 0
%                                   at the root to H.theta_tw at the tip)
%                    H.c         -> (data.nel x 1) chord at each beam node.
%                                   The cross-section at beam node i is
%                                   obtained by scaling the root section
%                                   by H.c(i)/H.c(1).
%   theta_0      : scalar root collective / pitch angle [rad]
%
%   Total section pitch about the span axis at beam node i is:
%       theta_total(i) = theta_0 + theta_tw(i) + theta_x(i)
%
%   The colorbar shows the vertical displacement of each skin vertex
%   relative to its undeformed (untapered-pose, but tapered-shape)
%   position:
%       w_v = z_deformed - z_undeformed
%   so the lift caused by torsional rotation is included.

    displacement = displacement(1:data.nel, :);
    theta_x      = 5*theta_x(1:data.nel);
    nnod = size(displacement, 1);   % number of beam (span) nodes

    % ---------------------------------------------------------------
    % 1. BUILD TOTAL PITCH AND GEOMETRIC TWIST PER BEAM NODE
    % ---------------------------------------------------------------
    theta_tw = H.theta_tw * linspace(0, 1, data.nel);

    % Pitch applied in the DEFORMED configuration (final orientation):
    theta_def = -theta_0 - theta_tw' - theta_x;        % nnod x 1

    % Pitch in the UNDEFORMED reference (no elastic torsion, no bending):
    theta_ref = -theta_0 - theta_tw';                  % nnod x 1

    % ---------------------------------------------------------------
    % 1b. CHORD-BASED TAPER SCALING PER BEAM NODE
    % ---------------------------------------------------------------
    % The provided cross-section mesh (nodes, elems) corresponds to the
    % ROOT element. The shape at any other span station is obtained by
    % scaling the root section by H.c(i)/H.c(1) (both y and z scale,
    % so the airfoil shrinks proportionally).
    c_vec = H.c(:);                       % column vector
    if numel(c_vec) ~= nnod
        % If H.c is given on a different grid, resample linearly
        c_vec = interp1(linspace(0, 1, numel(c_vec)), c_vec, ...
                        linspace(0, 1, nnod), 'linear')';
    end
    scale = c_vec / c_vec(1);             % nnod x 1   (scale(1) = 1)

    % ---------------------------------------------------------------
    % 2. EXTRACT ORDERED BOUNDARY LOOP OF THE CROSS-SECTION
    % ---------------------------------------------------------------
    all_edges = [];
    for e = 1:size(elems, 1)
        el = elems(e, :);
        el = el(el > 0);
        nv = numel(el);
        for k = 1:nv
            n1 = el(k);
            n2 = el(mod(k, nv) + 1);
            all_edges = [all_edges; sort([n1 n2])]; %#ok<AGROW>
        end
    end
    [uniq_edges, ~, ic] = unique(all_edges, 'rows');
    counts    = accumarray(ic, 1);
    bnd_edges = uniq_edges(counts == 1, :);

    % Walk the boundary into an ordered closed loop
    loop = zeros(size(bnd_edges, 1), 1);
    used = false(size(bnd_edges, 1), 1);
    loop(1) = bnd_edges(1, 1);
    loop(2) = bnd_edges(1, 2);
    used(1) = true;
    for k = 3:numel(loop)
        prev = loop(k - 1);
        idx  = find(~used & any(bnd_edges == prev, 2), 1);
        if isempty(idx)
            loop = loop(1:k-1);
            break
        end
        edge = bnd_edges(idx, :);
        loop(k) = edge(edge ~= prev);
        used(idx) = true;
    end
    nL = numel(loop);

    % ---------------------------------------------------------------
    % 3. UNDEFORMED REFERENCE BEAM-NODE POSITIONS
    % ---------------------------------------------------------------
    % We assume the beam axis is straight along X with the same X
    % coordinates as the deformed beam (axial stretch is negligible
    % for displacement-color purposes), and Y, Z equal to the root.
    x_ref = displacement(:, 1);
    y_ref = displacement(1, 2) * ones(nnod, 1); %#ok<NASGU>
    z_ref = displacement(1, 3) * ones(nnod, 1);

    % ---------------------------------------------------------------
    % 4. BUILD 3-D SKIN VERTICES (deformed) + UNDEFORMED Z FOR EACH
    % ---------------------------------------------------------------
    skin_xyz = zeros(nnod * nL, 3);
    skin_w   = zeros(nnod * nL, 1);

    for i = 1:nnod
        % Deformed pose
        x0 = displacement(i, 1);
        y0 = displacement(i, 2);
        z0 = displacement(i, 3);
        ct = cos(theta_def(i));
        st = sin(theta_def(i));

        % Reference (undeformed) pose
        z0r = z_ref(i);
        ctr = cos(theta_ref(i));
        str = sin(theta_ref(i));

        % Local taper scale at this beam node
        s = scale(i);

        for j = 1:nL
            cs_id = loop(j);
            % Scale the ROOT cross-section coordinates by the local
            % chord ratio to obtain the tapered local section.
            yl = s * nodes(cs_id, 1);
            zl = s * nodes(cs_id, 2);

            % Deformed vertex (rotation by theta_def about span axis,
            % then translated to deformed beam-node position)
            yr_d = ct * yl - st * zl;
            zr_d = st * yl + ct * zl;

            % Undeformed vertex Z (rotation by theta_ref only,
            % translated to root Z); the section is still tapered in
            % shape but sits at the undeformed beam position.
            zr_u = str * yl + ctr * zl;

            id = (i - 1) * nL + j;
            skin_xyz(id, 1) = x0;
            skin_xyz(id, 2) = y0 + yr_d;
            skin_xyz(id, 3) = z0 + zr_d;

            % Vertex vertical displacement = deformed Z - undeformed Z
            skin_w(id) = (z0 + zr_d) - (z0r + zr_u);
        end
    end

    % ---------------------------------------------------------------
    % 5. STITCH CONSECUTIVE CROSS-SECTIONS INTO QUAD STRIPS
    % ---------------------------------------------------------------
    nFaces = (nnod - 1) * nL;
    faces  = zeros(nFaces, 4);
    f = 1;
    for i = 1:(nnod - 1)
        o1 = (i - 1) * nL;
        o2 =  i      * nL;
        for j = 1:nL
            jn = mod(j, nL) + 1;
            faces(f, :) = [o1 + j, o1 + jn, o2 + jn, o2 + j];
            f = f + 1;
        end
    end

    % ---------------------------------------------------------------
    % 6. CAP THE ROOT AND TIP
    % ---------------------------------------------------------------
    root_face = (1:nL);
    tip_face  = (nnod - 1) * nL + (1:nL);

    % ---------------------------------------------------------------
    % 7. PLOT
    % ---------------------------------------------------------------
    blade = figure('Color', 'w');

    patch('Faces',           faces, ...
          'Vertices',        skin_xyz, ...
          'FaceVertexCData', skin_w, ...
          'FaceColor',       'interp', ...
          'EdgeColor',       'none');

    hold on
    patch('Faces',           root_face, ...
          'Vertices',        skin_xyz, ...
          'FaceVertexCData', skin_w, ...
          'FaceColor',       'interp', ...
          'EdgeColor',       'none');
    patch('Faces',           tip_face, ...
          'Vertices',        skin_xyz, ...
          'FaceVertexCData', skin_w, ...
          'FaceColor',       'interp', ...
          'EdgeColor',       'none');
    hold off

    axis equal
    view(3)
    grid on
    box on
    % Axes labels (LaTeX, size 18)
    xlabel('$X\,[m]$', 'Interpreter', 'latex', 'FontSize', 18)
    ylabel('$Y\,[m]$', 'Interpreter', 'latex', 'FontSize', 18)
    zlabel('$Z\,[m]$', 'Interpreter', 'latex', 'FontSize', 18)

    colormap(jet(256))
    cb = colorbar;
    ylabel(cb, 'Vertical displacement w [m]', 'Interpreter', 'latex', 'FontSize', 18)
    title('\textbf{Deformed Tapered Helicopter Blade - Vertical Displacement Distribution}', ...
          'Interpreter', 'latex', 'FontSize', 18)

    % Colorbar tick labels size
    cb.FontSize = 14;
end