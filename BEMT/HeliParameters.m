%% HELICOPTER DATA (Sikorsky UH-60 Black Hawk)
H = struct();

% Helicopter parameters
H.D = 16.36;                          % Rotor diameter [m]
H.A = H.D^2*pi/4;                     % Rotor surface [m^2]
H.W = 9800*9.8;                       % MTOW [N]

% Blade parameters
H.Nb    = 4;                          % Number of blades
H.Omega = 258;                        % Main rotor angular speed [rpm]
H.Vtip  = H.Omega*2*pi/60*H.D/2;      % Blade tip velocity [m/s]
H.e     = 0.75;                       % Oswald efficiency factor

H.Cl_alpha = 2*pi;                    % Lift slope
H.Cd0      = 0.008;                   % Parasite drag coefficient [-]
H.theta_tw = -18*2*pi/360;            % Blade geometrical twist [rad]
H.alph0    = 0;                       % Zero lift angle of attack [rad]

H.naca18   = 0.012;
H.naca10   = 0.007;

%% DENSITY DATA
Den = struct();

Den.h     = 0;                        % Height [m]
Den.rho0  = 1.225;                    % Density at 0 m [kg/m^3]
Den.T0    = 288.15;                   % Temperature at 0 m [K]
Den.R_air = 287.0528;                 % Specific gas constant [J⋅kg^−1⋅K^−1] 
Den.delta = 0.0065;                   % Lapse rate [K/m]