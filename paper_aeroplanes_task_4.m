%---- Paper Aeroplanes Task 4 ---------------------------------------------
%--------------------------------------------------------------------------

% close all previous plots
close all;
%----- Defining parameters ------------------------------------------------
D = [0, 1, sqrt(8), 4]; % parameter D
dt = 0.001; % setting time step
tspan = 0:dt:10; % time span
L = length(tspan);

%----- Setting up the figure with two subplots ----------------------------
figure('Units','normalized','OuterPosition',[0 0 1 1],'Theme','light');
figure('Units','normalized','OuterPosition',[0 0 1 1],'Theme','light');
% creating a matrix representing different colours for the plots 
% matrix has three columns representing the RGB values of the colour
colors = ['r', 'g', 'b', 'k'];

%----- Looping for each initial value of D --------------------------------
for i = 1:length(D)
    vel0 = [pi/4;
            1]; % initial conditions of (theta,s)'
%----- Defining the column vector that defines our system -----------------
    u = @(t,v) [(v(2)^2 - cos(v(1)))/v(2);
                -sin(v(1)) - D(i)*v(2)^2];

%----- Solving our system over the time domain tspan ----------------------
    % using Runge-Kutta 4 method
    vel = zeros(2,length(tspan)); vel(:,1) = [pi/4; 1];
    for m = 1:length(tspan)-1
        k1 = u(tspan(m),vel(:,m)); k2 = u(tspan(m),(vel(:,m) + 0.5*dt*k1));
        k3 = u(tspan(m),(vel(:,m) + 0.5*dt*k2)); k4 = u(tspan(m),(vel(:,m) + dt*k3));
        vel(:,m+1) = vel(:,m) + (dt/6)*(k1+2*k2+2*k3+k4);
    end

%----- Plotting the results -----------------------------------------------
    % adding to first subplot
    figure(1); hold on; grid on;
    % plotting theta over time
    plot(tspan, vel(1, :), '-', 'LineWidth', 2, 'Color', colors(i));
    % plotting speed over time
    plot(tspan, vel(2, :), '-.', 'LineWidth', 2, 'Color', colors(i));
    
    % colors(i,:) assigns one of the colours created in line 20 to a line
    % in our plot - want the lines to be the same colour for each iteration
    % so we know what lines correspond to each initial value of theta
    
%----- Plotting the trajectory on the xy-plane ----------------------------
    % adding to second subplot
    figure(2);
    hold on; grid on;
    % calculate x and y components of velocity through polar coordinates,
    % assuming initial position is at the origin
    dxdt = vel(2, :).*cos(vel(1, :));
    dydt = vel(2, :).*sin(vel(1, :));
    % finding x and y components of position by integrating
    x = cumtrapz(tspan, dxdt); % integrate to find x position
    y = cumtrapz(tspan, dydt); % integrate to find y position
    plot(x, y, '-', 'LineWidth', 2, 'Color', colors(i));

end

figure(1);
xlabel('Time (seconds)', 'Interpreter', 'latex');
ylabel('State Variables', 'Interpreter', 'latex');
legend('$\theta$,$D=0$', '$s$,$D=0$','$\theta$,$D=1$', '$s$,$D=1$','$\theta$,$D=\sqrt{8}$', '$s$,$D=\sqrt{8}$','$\theta$,$D=4$', '$s$,$D=4$', 'Interpreter', 'latex');
set(gca, 'FontSize', 44);
title('$\theta$ and $s$ over Time for $(\theta_0,s_0)=(\pi/4,1)$', 'Interpreter', 'latex');
ylim([-1.5 1.5]);

figure(2);
xlabel('X Position', 'Interpreter', 'latex');
ylabel('Y Position', 'Interpreter', 'latex');
title('Trajectories of Glider for different values of $D$','Interpreter','latex');
set(gca, 'FontSize', 44);
%legend('$D=0$', '$D=1$','$D=\sqrt{8}$', '$D=4$', 'Interpreter', 'latex');
ylim([-5.5 0.5]);
xlim([0 8.14]);