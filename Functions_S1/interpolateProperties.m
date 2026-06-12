function P = interpolateProperties(propTable, query_values, query_field)
% INTERPOLATEPROPERTIES  Linearly interpolate section properties from a
% pre-computed property table.
%
% INPUTS:
%   propTable    – struct array from generatePropertyTable (one entry per
%                  reference section). Must be sorted by the query field.
%   query_values – (N x 1) vector of values at which to interpolate
%                  (e.g., chord at each element midpoint, or NACA t/c).
%   query_field  – string: 'chord' or 'naca_tc' (the field used for
%                  interpolation abscissa).
%
% OUTPUT:
%   P – struct with (N x 1) vectors for each section property:
%       .A, .y_cg, .z_cg, .y_sc, .z_sc,
%       .Iy_cm, .Iz_cm, .Iyz_cm,
%       .Iy_c4, .Iz_c4, .Iyz_c4,
%       .J, .Isc, .ky, .kz, .kt
%
% The interpolation is LINEAR in the query variable. Values outside the
% table range are CLAMPED (no extrapolation).

N = length(query_values);

% Build the abscissa vector from the table
nRef  = length(propTable);
x_ref = zeros(nRef, 1);
for k = 1:nRef
    x_ref(k) = propTable(k).(query_field);
end

% Fields to interpolate (all scalar-valued per reference section)
fields = {'A', 'y_c4','y_cg', 'z_cg', 'y_sc', 'z_sc', ...
           'Iy_cm', 'Iz_cm', 'Iyz_cm', ...
           'Iy_sc', 'Iz_sc', 'Iyz_sc', ...
           'J', 'ky', 'kz', 'kt'};

% Clamp query values to table range
q = max(min(query_values, max(x_ref)), min(x_ref));

for f = 1:length(fields)
    fname = fields{f};

    % Build reference ordinate for this field
    y_ref = zeros(nRef, 1);
    for k = 1:nRef
        y_ref(k) = propTable(k).(fname);
    end

    % Interpolate
    P.(fname) = interp1(x_ref, y_ref, q, 'linear');
end

end
