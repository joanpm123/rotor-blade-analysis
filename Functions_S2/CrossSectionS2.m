function [sigma, tauy, tauz, VM] = CrossSectionS2(data, propTable_airfoil, secProp, ...
                                              ea, eb, es, et, r_break)
% CrossSectionS2  Internal stress recovery for Study 2 (airfoil variation).
%
% Each beam element uses either the root or tip cross-section node mesh
% depending on its spanwise position relative to r_break:
%   r <= r_break  →  root airfoil nodes
%   r >  r_break  →  tip  airfoil nodes
%
% INPUTS
%   data              – FEM data struct (field: nel)
%   propTable_airfoil – 1×2 struct; each entry has field .nodes [Ncs×2]
%   secProp           – struct with fields .E and .G (vectors of length nel)
%   ea, eb, es, et    – strain arrays from strainFunction
%   r_break           – non-dimensional spanwise break station ∈ (0,1]
%
% OUTPUTS
%   sigma  [Ncs_max × 1 × nel]  – normal stress  σxx  [Pa]
%   tauy   [Ncs_max × 1 × nel]  – shear stress   τxy  [Pa]
%   VM     [Ncs_max × 1 × nel]  – Von Mises      σVM  [Pa]

nodes_root = propTable_airfoil(1).nodes;   % [Ncs_root × 2]
nodes_tip  = propTable_airfoil(2).nodes;   % [Ncs_tip  × 2]

Ncs_max  = max(size(nodes_root, 1), size(nodes_tip, 1));
r_elem   = ((1:data.nel)' - 0.5) / data.nel;   % [nel × 1]

sigma = zeros(Ncs_max, 1, data.nel);
tauy  = zeros(Ncs_max, 1, data.nel);
VM    = zeros(Ncs_max, 1, data.nel);

for e = 1:data.nel

    % --- Select mesh based on radial position ---
    if r_elem(e) <= r_break
        nodes_e = nodes_root;
    else
        nodes_e = nodes_tip;
    end

    Ncs = size(nodes_e, 1);
    y_e = nodes_e(:, 1);
    z_e = nodes_e(:, 2);

    % --- Generalised strains ---
    exx = ea(1,e) + z_e * eb(1,e) - y_e * eb(2,e);
    gxy = es(1,e) - z_e * et(1,e);
    gxz = es(2,e) + y_e * et(1,e);

    % --- Material properties ---
    E_e = secProp.E(e);
    G_e = secProp.G(e);

    % --- Stresses ---
    sig_e = E_e * exx;
    txy_e = G_e * gxy;
    txz_e = G_e * gxz;

    % --- Von Mises ---
    VM_e = sqrt(sig_e.^2 + 3*(txy_e.^2 + txz_e.^2));

    % --- Store (zero-pad if this mesh is smaller than Ncs_max) ---
    sigma(1:Ncs, 1, e) = sig_e;
    tauy(1:Ncs,  1, e) = txy_e;
    tauz(1:Ncs,  1, e) = txz_e;
    VM(1:Ncs,    1, e) = VM_e;

end
end