%----- Paper Aeroplanes: A Dynamical Systems Perspective ------------------
%--------------------------------------------------------------------------
% close all previous plots
close all;
%----- Defining parameters ------------------------------------------------
D = 0; % parameter D
dt = 0.001; % setting time step
tspan = 0:dt:10; % time span
L = length(tspan);
% creating vector of different intial values of theta
theta_0 = [0, pi/8, pi/4];

%----- Defining the column vector that defines our system -----------------
u = @(t,v) [(v(2)^2 - cos(v(1)))/v(2);
            -sin(v(1)) - D*v(2)^2];

%----- Setting up two figures ---------------------------------------------
figure('Units','normalized','OuterPosition',[0 0 1 1],'Theme','light');
figure('Units','normalized','OuterPosition',[0 0 1 1],'Theme','light');
% creating a matrix representing different colours for the plots 
% matrix has three columns representing the RGB values of the colour
colors = ['r' 'g' 'b'];

%----- Looping for each initial value of theta ----------------------------
for i = 1:length(theta_0)
    vel0 = [theta_0(i);
            1]; % initial conditions of (theta,s)'

%----- Solving our system over the time domain tspan ----------------------
    opts = odeset('RelTol',1e-10,'AbsTol',1e-10);
    [t, vel] = ode45(@(t,v) u(t,v), tspan, vel0, opts);

%----- Plotting the results -----------------------------------------------
    % adding to first subplot
    figure(1); hold on; grid on;
    % plotting theta over time
    plot(t, vel(:, 1), '-', 'LineWidth', 2, 'Color', colors(i));
    % plotting speed over time
    plot(t, vel(:, 2), '-.', 'LineWidth', 2, 'Color', colors(i));
    
    % colors(i,:) assigns one of the colours created in line 20 to a line
    % in our plot - want the lines to be the same colour for each iteration
    % so we know what lines correspond to each initial value of theta
    
%----- Plotting the trajectory on the xy-plane ----------------------------
    % adding to second subplot
    figure(2);
    hold on; grid on;
    % calculate x and y components of velocity through polar coordinates,
    % assuming initial position is at the origin
    dxdt = vel(:, 2).*cos(vel(:, 1));
    dydt = vel(:, 2).*sin(vel(:, 1));
    % finding x and y components of position by integrating
    x = cumtrapz(t, dxdt); % integrate to find x position
    y = cumtrapz(t, dydt); % integrate to find y position
    plot(x, y, '-', 'LineWidth', 2, 'Color', colors(i));
    xlabel('X Position', 'Interpreter', 'latex');
    ylabel('Y Position', 'Interpreter', 'latex');
    %legend('show', 'Interpreter', 'latex');
    title('Trajectories of Glider for different values of $\theta_0$','Interpreter','latex');
    set(gca, 'FontSize', 48);

end

figure(1);
% creating labels
ylim([-1 1.5]);
xlabel('Time (seconds)', 'Interpreter', 'latex');
ylabel('State Variables', 'Interpreter', 'latex');
legend('$\theta, \theta_0=0$', '$s, \theta_0=0$', '$\theta, \theta_0=\pi/8$', '$s, \theta_0=\pi/8$', '$\theta, \theta_0=\pi/4$', '$s, \theta_0=\pi/4$', 'Interpreter', 'latex');
set(gca, 'FontSize', 48);
title('$\theta$ and $s$ over Time', 'Interpreter', 'latex');

figure(2);
%legend('$\theta, \theta(0)=0$','$\theta, \theta(0)=\pi/8$','$\theta, \theta(0)=\pi/4$', 'Interpreter', 'latex');
ylim([-0.7 0.5]);
