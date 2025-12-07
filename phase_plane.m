close all;
%set D for plot
%PLOTS FOR D 1 TO 4, EQUILIBRIUM POINTS
D = [0, 1, 2, 3, 4];
%creates coordinate grid
[X, Y] = meshgrid(-pi/2:pi/20:pi/2,0:0.05:2);
%plot breaks with scaling when step divides domain
%scale factor for plotting
sf = 0.1;
%equilibrium points for each value of D, calculated analytically beforehand
equilibrium_theta = [atan(0),atan(-1),atan(-2),atan(-3),atan(-4)];
equilibrium_s = sqrt(cos(equilibrium_theta));
for i = 1:length(D)
    figure('Theme', 'light');
    %set up functions from X,Y -> U,V
    U = (Y.^2 - cos(X))./Y; %normalise vectors to length 1, scaling afterwards. show length with colours?
    V = -sin(X) - D(i)*Y.^2;
    %setup grids for U direction, V direction, and Magnitude
    Unorm = zeros(length(X),length(Y));
    Vnorm = zeros(length(X),length(Y));
    magnitude = hypot(U,V);
    %Scale U and V
    Unorm = sf.*U./magnitude;
    Vnorm = sf.*V./magnitude;
    hold on;
    %Plot
    quiver(X,Y,Unorm,Vnorm, 'HandleVisibility', 'off', 'Color', 'black');
    %shows how system develops from initial conditions
    streamline(X,Y,U,V,0,1, 'LineWidth', 2, 'DisplayName', '$(s_0,\theta_0)=(1,0)$', 'Color', 'r');
    streamline(X,Y,U,V,pi/8,1, 'LineWidth', 2, 'DisplayName', '$(s_0,\theta_0)=(1,\pi/8)$', 'Color', 'g');
    streamline(X,Y,U,V,pi/4,1, 'LineWidth', 2, 'DisplayName', '$(s_0,\theta_0)=(1,\pi/4)$', 'Color', 'b');
    %PLOT EQUILIBRIUM POINTS MANUALLY LIKE IN TASK 3
    plot(equilibrium_theta(i), equilibrium_s(i), '.', 'Color', 'm', 'LineWidth', 2, 'HandleVisibility', 'off','DisplayName','Equilibrium Point');
    %axis equal;
    xlabel('$\theta$', 'Interpreter', 'latex','Position',[1.65,0.085,0]);
    ylabel('$s$', 'Interpreter', 'latex', 'Rotation', 0,'Position',[-1.57,1.1,0]);
    set(gca, 'FontSize', 48);
    axis([-pi/2 pi/2 0 2]);
    %legend('show','Location', 'best', 'Interpreter', 'latex');
    %title(sprintf('Phase Plane ($D = %i$)', i-1), 'Interpreter', 'latex');
end

figure(1);
xlabel('$\theta$', 'Interpreter', 'latex','Position',[1.65,0.085,0]);
ylabel('$s$', 'Interpreter', 'latex', 'Rotation', 0,'Position',[-1.57,1.75,0]);
ylim([0 1.75]);
figure(2);
legend('show','Location', 'best', 'Interpreter', 'latex');
ylim([0 1.1]);
figure(3);
ylim([0 1.1]);
figure(4);
legend('show','Location', 'best', 'Interpreter', 'latex');
ylim([0 1.1]);
figure(5);
ylim([0 1.1]);