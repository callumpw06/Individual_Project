% Task 3
D_values = linspace(0,4,400);
s_star   = (1 + D_values.^2).^(-1/4);%equilibrium speed
D_star = sqrt(8);%bifurcation point
theta_star = -atan(D_values);

%Each different region           
spiral = (D_values > 0) & (D_values < D_star);
node   = (D_values > D_star);

%Plotting Values of S
figure('Theme','light'); hold on;
%centre
plot(0,1,'ko','MarkerFaceColor','k','DisplayName','Bifurcation point ($D = 0$)','MarkerSize',14);
%spiral sink
plot(D_values(spiral),s_star(spiral),'b','LineWidth',1.8,'DisplayName','Spiral sink ($0 < D < \sqrt{8}$)');
%nodal sink
plot(D_values(node),s_star(node),'r','LineWidth',1.8,'DisplayName','Nodal sink ($D > \sqrt{8}$)');
% mark the bifurcation point as D = sqrt(8)
s_star = (1+ D_star^2)^(-1/4);
plot(D_star, s_star, 'ms', 'MarkerFaceColor','m','DisplayName','Bifurcation point (D = $\sqrt{8}$)','MarkerSize',14);

xlabel('$D$', 'Interpreter', 'latex');
ylabel('Equilibrium speed $s^*$', 'Interpreter', 'latex');
title('Bifurcation diagram of $D$ against $s^*$', 'Interpreter','latex');
grid on;
%legend('Location','best','Interpreter','latex');
set(gca, 'FontSize', 44);

%Plotting Values of Theta
figure('Theme','light'); hold on;
%centre
plot(0,0,'ko','MarkerFaceColor','k','DisplayName','Bifurcation point ($D = 0$)','MarkerSize',14);
%spiral sink
plot(D_values(spiral),theta_star(spiral),'b','LineWidth',1.8,'DisplayName','Spiral sink ($0 < D < \sqrt{8}$)');
%nodal sink
plot(D_values(node),theta_star(node),'r','LineWidth',1.8,'DisplayName','Nodal sink ($D > \sqrt{8}$)');
% mark the bifurcation point as D = sqrt(8)
theta_star = -atan(D_star);
plot(D_star, theta_star, 'ms', 'MarkerFaceColor','m','DisplayName','Bifurcation point (D = $\sqrt{8}$)','MarkerSize',14);

xlabel('$D$', 'Interpreter', 'latex');
ylabel('Equilibrium AoA $\theta^*$', 'Interpreter', 'latex');
title('Bifurcation diagram of $D$ against $\theta^*$', 'Interpreter','latex');
grid on;
legend('Location','best','Interpreter','latex');
set(gca, 'FontSize', 44);
