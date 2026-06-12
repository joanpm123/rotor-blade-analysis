function P = interpAirfoil(propTable_airfoil, r_elem, r_break)
% INTERPAIRFOIL  Linear interpolation between root and tip airfoil properties.
%   r < r_break*0.8 : pure root airfoil
%   r > r_break*1.2 : pure tip airfoil
%   Between: linear transition zone centred at r_break.

    % Blending factor: 0 at root, 1 at tip, transition around r_break
    alpha = max(0, min(1, (r_elem - r_break*0.8) ./ (r_break*1.2 - r_break*0.8 + 1e-10)));

    fields = {'A', 'y_cg', 'z_cg', 'y_c4', 'z_c4', ...
              'Iy_cm', 'Iz_cm', 'Iyz_cm', ...
              'Iy_sc', 'Iz_sc', 'Iyz_sc', ...
              'J', 'ky', 'kz', 'kt'};

    P = struct();
    for f = 1:length(fields)
        fname = fields{f};
        val_root = propTable_airfoil(1).(fname);
        val_tip  = propTable_airfoil(2).(fname);
        P.(fname) = val_root * (1 - alpha) + val_tip * alpha;
    end
end