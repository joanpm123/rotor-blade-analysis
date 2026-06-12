function [Kel,Ka,Ba,Kb,Bb,Ks,Bs,Kt,Bt,R,T] = stiffnessFunction(data,x,Tn,secProp,j)

% Matrices initiation
Kel = zeros(data.nne*data.ni,data.nne*data.ni,data.nel);
Ka  = zeros(data.nne*data.ni,data.nne*data.ni,data.nel);
Kb  = zeros(data.nne*data.ni,data.nne*data.ni,data.nel);
Ks  = zeros(data.nne*data.ni,data.nne*data.ni,data.nel);
Kt  = zeros(data.nne*data.ni,data.nne*data.ni,data.nel);
R  = zeros(data.nne*data.ni,data.nne*data.ni,data.nel);
Ba = cell(data.nel,1);
Bb = cell(data.nel,1);
Bs = cell(data.nel,1);
Bt = cell(data.nel,1);

% Loop over elements
for e = 1:data.nel
    l = norm(x(Tn(e,2),:)-x(Tn(e,1),:));     % Element size

    % Element local axis
    i = (transpose(x(Tn(e,2),:))-transpose(x(Tn(e,1),:)))/l; % Exp.(4.74)
    k = cross(i,j).';

    % Rotation matrix
    z = zeros(3,3);
    k = k';
    r = [i j k z;
         z i j k]';
    Z = zeros(6,6);
    R(:,:,e) = [r Z;
                Z r];
    
    % Centres
    y_cg = secProp.y_cg(e);
    z_cg = secProp.z_cg(e);
    y_sc = 0;
    z_sc = 0;

    % Translation matrix
    t = [1, 0, 0,     0, z_cg, -y_cg;
         0, 1, 0, -z_sc,    0,     0;
         0, 0, 1,  y_sc,    0,     0;
         0, 0, 0,     1,    0,     0;
         0, 0, 0,     0,    1,     0;
         0, 0, 0,     0,    0,     1];
    
    T(:,:,e) =[t Z;
               Z t];


    % Element shape function derivatives
    dN(1) = -1/l;
    dN(2) = 1/l;

    % Axial component
    [Ka(:,:,e),Ba{e}] = axialFunction(data,l,R(:,:,e),dN,secProp,e);

    % Bending component
    [Kb(:,:,e),Bb{e}] = bendingFunction(data,l,R(:,:,e),dN,secProp,e);

    % Shear component
    [Ks(:,:,e),Bs{e}] = shearFunction(data,l,R(:,:,e),dN,secProp,e);

    % Torsion component
    [Kt(:,:,e),Bt{e}] = torsionFunction(data,l,R(:,:,e),dN,secProp,e);

    % Element stiffness matrix
    %Kel(:,:,e) = Ka(:,:,e) + Kb(:,:,e) + Ks(:,:,e) + Kt(:,:,e);

    % Constitutive matrix
    E   = secProp.E(e);
    G   = secProp.G(e);
    A   = secProp.A(e);
    ky  = secProp.ky(e);
    kz  = secProp.kz(e);
    kt  = secProp.kt(e);
    Iy  = secProp.Iy(e);
    Iz  = secProp.Iz(e);
    Iyz = secProp.Iyz(e);
    J   = secProp.J(e);
    
    %
    C = [ E*A,         0,        0,      0,     0,     0;
            0,    ky*G*A,        0,      0,     0,     0;
            0,         0,   kz*G*A,      0,     0,     0;
            0,         0,        0, kt*G*J,     0,     0;
            0,         0,        0,      0,  E*Iy, E*Iyz;
            0,         0,        0,      0, E*Iyz,  E*Iz];

    B = [Ba{e}; Bs{e}; Bt{e}; Bb{e}];

    Kel(:,:,e) = l*R(:,:,e)'*T(:,:,e)*B'*C*B*T(:,:,e)*R(:,:,e);
    
end

end