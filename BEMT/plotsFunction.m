function fig2 = plotsFunction(rj,dCT,lambda_c)

% Plotting dC_T/dr distributions
fig2 = figure;
box on 
grid on 
title('dC_T/dr distributions')
ylabel('dC_T/dr')
xlabel('r')
hold on 
plot(rj,dCT(:,1),'LineWidth', 2 , 'DisplayName',sprintf('\\lambda_c= %.2f',lambda_c))
legend('Location','Best') 
hold off
end