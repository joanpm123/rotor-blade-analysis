function blade = BladePlotTaper(displacement, sigma, elems, nodes, data, t, theta_x, H, theta_0)
% BLADEPLOTTAPER  Render the deformed TAPERED helicopter blade with
%                 per-vertex stress.
%
%   displacement : (data.nel x 3)              [x y z] of each beam node
%   sigma        : (size(nodes,1) x 1 x data.nel)  stress at every
%                                                  cross-section node, for
%                                                  every beam node
%   elems        : (ne x k)   ROOT cross-section connectivity (zero-padded)
%   nodes        : (ncs x 2)  ROOT cross-section [y z] coordinates
%   data.nel     : number of beam nodes
%   t            : 1 -> sigma_xx title, otherwise -> Von Mises title
%   theta_x      : (data.nel x 1) elastic torsion angle [rad] at each beam node
%   H            : struct containing
%                    H.theta_tw  -> tip geometric twist [rad] (scalar; the
%                                   twist is distributed linearly from 0
%                                   at the root to H.theta_tw at the tip)
%                    H.c         -> (data.nel x 1) chord at each beam node.
%                                   The cross-section at beam node i is
%                                   obtained by scaling the root section
%                                   by H.c(i)/H.c(1).
%   theta_0      : scalar root collective / pitch angle [rad]

    displacement = displacement(1:data.nel, :);
    theta_x      = 1*(theta_x(1:data.nel));
    nnod = size(displacement, 1);   % number of beam (span) nodes
    ncs  = size(nodes, 1);          % number of cross-section nodes %#ok<NASGU>


    % ---------------------------------------------------------------
    % 0. BUILD TOTAL PITCH AND GEOMETRIC TWIST PER BEAM NODE
    % ---------------------------------------------------------------
    theta_tw = H.theta_tw * linspace(0, 1, data.nel);

    % Pitch applied in the DEFORMED configuration (final orientation):
    theta_def = -theta_0 - theta_tw' + theta_x;        % nnod x 1

    % Pitch in the UNDEFORMED reference (no elastic torsion, no bending):
    theta_ref = -theta_0 - theta_tw';                  % nnod x 1


    % ---------------------------------------------------------------
    % 0b. CHORD-BASED TAPER SCALING PER BEAM NODE
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
    % 1. EXTRACT ORDERED BOUNDARY LOOP OF THE CROSS-SECTION
    % ---------------------------------------------------------------
    % First find the unordered boundary edges (edges that belong to a
    % single 2-D element), then walk them to build a closed loop. The
    % ordered loop lets us stitch consecutive cross-sections together
    % into clean quad strips with no degenerate or crossed faces.

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
    counts = accumarray(ic, 1);
    bnd_edges = uniq_edges(counts == 1, :);   % each row is one edge

    % Walk the boundary to obtain an ordered list of node indices
    loop = zeros(size(bnd_edges, 1), 1);
    used = false(size(bnd_edges, 1), 1);
    loop(1) = bnd_edges(1, 1);
    loop(2) = bnd_edges(1, 2);
    used(1) = true;
    for k = 3:numel(loop)
        prev = loop(k - 1);
        % find the unused edge that contains 'prev'
        idx = find(~used & any(bnd_edges == prev, 2), 1);
        if isempty(idx)
            % open boundary (shouldn't happen for closed airfoil)
            loop = loop(1:k-1);
            break
        end
        edge = bnd_edges(idx, :);
        loop(k) = edge(edge ~= prev);
        used(idx) = true;
    end
    nL = numel(loop);

    % We assume the beam axis is straight along X with the same X
    % coordinates as the deformed beam (axial stretch is negligible
    % for displacement-color purposes), and Y, Z equal to the root.
    x_ref = displacement(:, 1);                       %#ok<NASGU>
    y_ref = displacement(1, 2) * ones(nnod, 1);       %#ok<NASGU>
    z_ref = displacement(1, 3) * ones(nnod, 1);       %#ok<NASGU>


    % ---------------------------------------------------------------
    % 2. BUILD 3-D SKIN VERTICES (only boundary nodes are needed)
    % ---------------------------------------------------------------
    skin_xyz   = zeros(nnod * nL, 3);
    skin_sigma = zeros(nnod * nL, 1);

    for i = 1:nnod
        % Deformed pose
        x0 = displacement(i, 1);
        y0 = displacement(i, 2);
        z0 = displacement(i, 3);
        ct = cos(theta_def(i));
        st = sin(theta_def(i));

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

            id = (i - 1) * nL + j;
            skin_xyz(id, 1) = x0;
            skin_xyz(id, 2) = y0 + yr_d;
            skin_xyz(id, 3) = z0 + zr_d;
            skin_sigma(id)  = sigma(cs_id, 1, i);
        end
    end


    % ---------------------------------------------------------------
    % 3. STITCH CONSECUTIVE CROSS-SECTIONS INTO QUAD STRIPS
    % ---------------------------------------------------------------
    nFaces = (nnod - 1) * nL;
    faces  = zeros(nFaces, 4);
    f = 1;
    for i = 1:(nnod - 1)
        o1 = (i - 1) * nL;
        o2 =  i      * nL;
        for j = 1:nL
            jn = mod(j, nL) + 1;          % wrap around the loop
            faces(f, :) = [o1 + j, o1 + jn, o2 + jn, o2 + j];
            f = f + 1;
        end
    end

    % ---------------------------------------------------------------
    % 4. CAP THE ROOT AND THE TIP
    % ---------------------------------------------------------------
    % Use the ordered loop as a single n-gon face so the blade looks
    % closed at both ends.
    root_face = (1:nL);                   % first cross-section
    tip_face  = (nnod - 1) * nL + (1:nL); % last cross-section

    % ---------------------------------------------------------------
    % 5. PLOT
    % ---------------------------------------------------------------
    blade = figure('Color', 'w');

    % Skin (smooth per-vertex stress)
    patch('Faces',          faces, ...
          'Vertices',       skin_xyz, ...
          'FaceVertexCData', skin_sigma, ...
          'FaceColor',     'interp', ...
          'EdgeColor',     'none');

    hold on

    % End caps - same colormap, no edges
    patch('Faces',          root_face, ...
          'Vertices',       skin_xyz, ...
          'FaceVertexCData', skin_sigma, ...
          'FaceColor',     'interp', ...
          'EdgeColor',     'none');
    patch('Faces',          tip_face, ...
          'Vertices',       skin_xyz, ...
          'FaceVertexCData', skin_sigma, ...
          'FaceColor',     'interp', ...
          'EdgeColor',     'none');

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

    % Colorbar label (LaTeX)
    if t == 1
        ylabel(cb, '$\sigma_{xx}\,[Pa]$', 'Interpreter', 'latex', 'FontSize', 18)
        title('\textbf{Deformed Tapered Helicopter Blade - $\sigma_{xx}$ Distribution}', ...
              'Interpreter', 'latex', 'FontSize', 18)
    else
        ylabel(cb, '$\sigma_{VM}\,[Pa]$', 'Interpreter', 'latex', 'FontSize', 18)
        title('\textbf{Deformed Tapered Helicopter Blade - Von Mises Stress Distribution}', ...
              'Interpreter', 'latex', 'FontSize', 18)
    end

    % Colorbar tick labels size
    cb.FontSize = 14;

end