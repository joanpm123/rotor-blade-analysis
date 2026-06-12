%% 3D Cantilever Beam Geometry (Timoshenko Beam)

% Helicopter and Blade parameters
HeliParameters
% -----------------------
% Parameters
% -----------------------
L  = H.D/2;        % Total beam length [m]
ne = 150;            % Number of elements
nn = ne + 1;       % Number of nodes


% Beam direction (arbitrary 3D orientation)
d = [2, 0, 0];
d = d / norm(d);     % Unit direction vector

% Node spacing
ds = L / ne;

% -----------------------
% Nodal coordinates
% -----------------------
coordinates = zeros(nn,3);

for i = 1:nn
    coordinates(i,:) = (i-1) * ds * d;
end

% -----------------------
% Element connectivity
% -----------------------
% Each row: [node_i  node_j]
connectivity = zeros(ne,2);

for e = 1:ne
    connectivity(e,:) = [e, e+1];
end

% Display
disp('Node coordinates (x y z):')
disp(coordinates)

disp('Element connectivity:')
disp(connectivity)
