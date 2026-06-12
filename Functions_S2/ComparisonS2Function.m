function [fig1,fig2,fig3,fig4,fig5,fig6] = ComparisonS2Function(results,nCases)

%------------------------------------------------------------------------%

fig3 = figure('Name','Integrated Induced and Profile Powers - Chord Study');

hold on;
grid on;

% Preallocate
Pi = zeros(1,nCases);
P0 = zeros(1,nCases);
rbreak = zeros(1,nCases);

for iCase = 1:nCases  
    % Integrated induced power
    Pi(iCase) = sum(results(iCase).dCP_i);

    % Integrated profile power
    P0(iCase) = sum(results(iCase).dCP_0);

    % rbreak ratio
    rbreak(iCase) = results(iCase).r_break;
end

% Optional: sort by rbreak ratio
[rbreak, idx] = sort(rbreak);

Pi = Pi(idx);
P0 = P0(idx);

%% LEFT AXIS : INDUCED POWER
yyaxis left

plot(rbreak, Pi, '-o','LineWidth', 2,'MarkerSize', 7);
set(gca, 'FontSize', 12);
ylabel('$P_i$ [W]','Interpreter', 'latex','FontSize', 16.5);

%% RIGHT AXIS : PROFILE POWER
yyaxis right

plot(rbreak, P0, '-s', ...
    'LineWidth', 2, ...
    'MarkerSize', 7);

ylabel('$P_0$ [W]','Interpreter', 'latex','FontSize', 16.5);

%% COMMON FORMATTING
xlabel('Transition point [$r_b$]', 'Interpreter', 'latex', 'FontSize', 16.5);

title('\textbf{Power Variation with Airfoil Transition}', 'Interpreter', 'latex', ...
    'FontSize', 18);

legend({'$P_i$', '$P_0$'}, 'Interpreter', 'latex', 'FontSize', 16, 'Location', 'best');

hold off;

%-------------------------------------------------------------------------%
fig4 = figure('Name','Figure of Merit - Chord Study');

hold on;
grid on;

FM = zeros(1,nCases);
for iCase = 1:nCases  
    FM(iCase) = results(iCase).FM;
end

plot(rbreak, FM, '-o','LineWidth', 2,'MarkerSize', 7);

set(gca, 'FontSize', 12);
title('$\textbf{Figure of Merit Variation with Airfoil Transition}$', 'Interpreter', 'latex', 'FontSize', 18)
xlabel('Transition point $[r_b]$', 'Interpreter', 'latex', 'FontSize', 16.5)
ylabel('FM', 'Interpreter', 'latex', 'FontSize', 16.5)



%-------------------------------------------------------------------------%
w_tip = zeros(1,nCases);
twist = zeros(1,nCases);
rot_y = zeros(1,nCases);

for iCase = 1:nCases
    w_tip(iCase) = results(iCase).w_tip;
    twist(iCase) = rad2deg(results(iCase).twist_tip);
    rot_y(iCase) = abs(max(rad2deg(results(iCase).theta_y)));
end


fig2 = figure('Name', 'Flapwise Deflection - Chord Study');

hold on;
plot(rbreak,w_tip,'-o','LineWidth', 2,'MarkerSize', 7);
xlim([min(rbreak) max(rbreak)])
grid on;
set(gca, 'FontSize', 12);

title('$\textbf{Blade Tip Vertical Displacement}$', 'Interpreter', 'latex', 'FontSize', 18)
xlabel('Transition point [$r_b$]', 'Interpreter', 'latex', 'FontSize', 16.5)
ylabel('Deflection [m]', 'Interpreter', 'latex', 'FontSize', 16.5)

%------------------------------------------------------------------------%
fig5 = figure('Name', 'Angles - Chord Study');
hold on;
plot(rbreak,twist,'LineWidth', 2,'MarkerSize', 7);
plot(rbreak,rot_y,'LineWidth', 2,'MarkerSize', 7);
xlim([min(rbreak) max(rbreak)])
grid on;
set(gca, 'FontSize', 12);

title('$\textbf{Maximum Rotation Angles}$', 'Interpreter', 'latex', 'FontSize', 18)
xlabel('Transition point [$r_b$]', 'Interpreter', 'latex', 'FontSize', 16.5)
ylabel('Angle [deg]', 'Interpreter', 'latex', 'FontSize', 16.5)
legend({'$\theta_x$', '$\theta_y$'}, 'Interpreter', 'latex', 'FontSize', 16, 'Location', 'best');

%------------------------------------------------------------------------%
Mx = zeros(1,nCases);
My = zeros(1,nCases);
Nx = zeros(1,nCases);

for iCase = 1:nCases
    Mx(iCase) = -results(iCase).Mx(1)/1000;
    My(iCase) = results(iCase).My(1)/1000;
    Nx(iCase) = results(iCase).Fx(1)/1000;
end

fig1 = figure('Name', 'Moments - Chord Study');
hold on;
plot(rbreak,Mx,'LineWidth', 2,'MarkerSize', 7);
plot(rbreak,My,'LineWidth', 2,'MarkerSize', 7);
xlim([min(rbreak) max(rbreak)])
grid on;
set(gca, 'FontSize', 12);

title('$\textbf{Root Flapwise Bending and Torsional Moments}$', 'Interpreter', 'latex', 'FontSize', 18)
xlabel('Transition point [$r_b$]', 'Interpreter', 'latex', 'FontSize', 16.5)
ylabel('Moment [kN m]', 'Interpreter', 'latex', 'FontSize', 16.5)
legend({'$M_x$', '$M_y$'}, 'Interpreter', 'latex', 'FontSize', 16, 'Location', 'best');

%------------------------------------------------------------------------%
fig6 = figure('Name', 'Axial Force - Chord Study');
hold on;
plot(rbreak,Nx,'LineWidth', 2,'MarkerSize', 7);
xlim([min(rbreak) max(rbreak)])
grid on;
set(gca, 'FontSize', 12);

title('$\textbf{Root Axial Force Variation}$', 'Interpreter', 'latex', 'FontSize', 18)
xlabel('Transition point [$r_b$]', 'Interpreter', 'latex', 'FontSize', 16.5)
ylabel('Force [kN]', 'Interpreter', 'latex', 'FontSize', 16.5)

end