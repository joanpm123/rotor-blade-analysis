function [figure1,figure2,figure3,figure4] = figuresFunction(x,ue,displacement,Fx,Fy,Fz,Mx,My,Mz)

figure1 = [];

%-------------------------------------------------%
displacement = [ue(1:3,:)' ; ue(7:9,end)'];
figure2 = figure;
subplot(3,1,1)
plot(x(:,1),displacement(:,1),'b','LineWidth',2)
xlim([x(1,1), x(end,1)]);
set(gca, 'FontSize', 10);
title('$\textbf{Displacement X}$', 'Interpreter', 'latex', 'FontSize', 11)
ylabel('Amplitude', 'Interpreter', 'latex', 'FontSize', 11)
grid on
subplot(3,1,2)
plot(x(:,1),displacement(:,2),'r','LineWidth',2)
xlim([x(1,1), x(end,1)]);
set(gca, 'FontSize', 10);
title('$\textbf{Displacement Y}$', 'Interpreter', 'latex', 'FontSize', 11)
ylabel('Amplitude', 'Interpreter', 'latex', 'FontSize', 11)
grid on
subplot(3,1,3)
plot(x(:,1),displacement(:,3),'k','LineWidth',2)
xlim([x(1,1), x(end,1)]);
set(gca, 'FontSize', 10);
title('$\textbf{Displacement Z}$', 'Interpreter', 'latex', 'FontSize', 11)
xlabel('Index', 'Interpreter', 'latex', 'FontSize', 11)
ylabel('Amplitude', 'Interpreter', 'latex', 'FontSize', 11)
grid on
%-------------------------------------------------%
figure3 = figure;
subplot(3,1,1)
plot(x(1:end-1,1),Fx,'b','LineWidth',2)
xlim([x(1,1), x(end,1)]);
set(gca, 'FontSize', 10);
title('\textbf{Axial force $F_x$}', 'Interpreter', 'latex', 'FontSize', 11)
xlabel('Span x [m]', 'Interpreter', 'latex', 'FontSize', 11)
ylabel('Force [N]', 'Interpreter', 'latex', 'FontSize', 11)
grid on

%-------------------------------------------------%
figure4 = figure;
subplot(3,1,1)
plot(x(1:end-1,1),-Mx,'b','LineWidth',2)
xlim([x(1,1), x(end,1)]);
set(gca, 'FontSize', 10);
title('\textbf{Torsional moment $M_x$}', 'Interpreter', 'latex', 'FontSize', 11)
ylabel('Moment [Nm]', 'Interpreter', 'latex', 'FontSize', 11)
grid on
subplot(3,1,2)
plot(x(1:end-1,1),My,'r','LineWidth',2)
xlim([x(1,1), x(end,1)]);
set(gca, 'FontSize', 10);
title('\textbf{Bending moment $M_y$}', 'Interpreter', 'latex', 'FontSize', 11)
ylabel('Moment [Nm]', 'Interpreter', 'latex', 'FontSize', 11)
grid on
subplot(3,1,3)
plot(x(1:end-1,1),Mz,'k','LineWidth',2)
xlim([x(1,1), x(end,1)]);
set(gca, 'FontSize', 10);
title('\textbf{Lagging moment $M_z$}', 'Interpreter', 'latex', 'FontSize', 11)
xlabel('Span x [m]', 'Interpreter', 'latex', 'FontSize', 11)
ylabel('Moment [Nm]', 'Interpreter', 'latex', 'FontSize', 11)
grid on
end