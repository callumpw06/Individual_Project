%----- Paper Aeroplanes: A Dynamical Systems Perspective ------------------
%--------------------------------------------------------------------------
% close all previous plots
close all;
%----- Defining parameters ------------------------------------------------
D = 0; % parameter D
dt = 0.001; % setting time step
tspan = 0:dt:10; % time span
L = length(tspan);
% assigning the intial value of theta
theta_0 = 0;
%----- Defining the column vector that defines our system -----------------
%different u's for different values of D
u0 = @(t,v) [(v(2)^2 - cos(v(1)))/v(2);
            -sin(v(1)) - 0*v(2)^2];
u1 = @(t,v) [(v(2)^2 - cos(v(1)))/v(2);
            -sin(v(1)) - 1*v(2)^2];
u2 = @(t,v) [(v(2)^2 - cos(v(1)))/v(2);
            -sin(v(1)) - 2*v(2)^2];
u3 = @(t,v) [(v(2)^2 - cos(v(1)))/v(2);
            -sin(v(1)) - 3*v(2)^2];
u4 = @(t,v) [(v(2)^2 - cos(v(1)))/v(2);
            -sin(v(1)) - 4*v(2)^2];
%----- Setting up the figure with two subplots ----------------------------
figure('Units','normalized','OuterPosition',[0 0 1 1],'Theme','light');
figure('Units','normalized','OuterPosition',[0 0 1 1],'Theme','light');
% creating a matrix representing different colours for the plots 
% matrix has three columns representing the RGB values of the colour

%----- Assign conditions for the initial value of theta ----------------------------
    init_vel0 = [theta_0; 1];
% initial conditions of (theta,s)'

%----- Solving our systems over the time domain tspan ----------------------
    opts = odeset('RelTol',1e-10,'AbsTol',1e-10);
    [t, vel0] = ode45(@(t,v) u0(t,v), tspan, init_vel0, opts);
    [t, vel1] = ode45(@(t,v) u1(t,v), tspan, init_vel0, opts);
    [t, vel2] = ode45(@(t,v) u2(t,v), tspan, init_vel0, opts);
    [t, vel3] = ode45(@(t,v) u3(t,v), tspan, init_vel0, opts);
    [t, vel4] = ode45(@(t,v) u4(t,v), tspan, init_vel0, opts);

%----- Plotting the results -----------------------------------------------
    % adding to first subplot
    figure(1); hold on; grid on;
    % plotting theta over time
    plot(t, vel0(:, 1), '-r', 'LineWidth', 2);
    % plotting speed over time
    plot(t, vel0(:, 2), '-.r', 'LineWidth', 2);
    % plotting theta over time
    plot(t, vel1(:, 1), '-g', 'LineWidth', 2);
    % plotting speed over time
    plot(t, vel1(:, 2), '-.g', 'LineWidth', 2);
    % plotting theta over time
    plot(t, vel2(:, 1), '-b', 'LineWidth', 2);
    % plotting speed over time
    plot(t, vel2(:, 2), '-.b', 'LineWidth', 2);
    % plotting theta over time
    plot(t, vel3(:, 1), '-m', 'LineWidth', 2);
    % plotting speed over time
    plot(t, vel3(:, 2), '-.m', 'LineWidth', 2);
    % plotting theta over time
    plot(t, vel4(:, 1), '-k', 'LineWidth', 2);
    % plotting speed over time
    plot(t, vel4(:, 2), '-.k', 'LineWidth', 2);
    % creating labels
xlabel('Time (seconds)', 'Interpreter', 'latex');
ylabel('State Variables', 'Interpreter', 'latex');
legend('$\theta$,($D = 0$)','$s$,($D = 0$)','$\theta$,($D = 1$)','$s$,($D = 1$)','$\theta$,($D = 2$)','$s$,($D = 2$)','$\theta$,($D = 3$)','$s$,($D = 3$)','$\theta$,($D = 4$)','$s$,($D = 4$)','Interpreter', 'latex');
set(gca, 'FontSize', 40);

title('Dynamics of Paper Aeroplanes', 'Interpreter', 'latex');
    
%----- Plotting the trajectory on the xy-plane ----------------------------
    % adding to second subplot
    figure(2);
    hold on; grid on;
    % calculate x and y components of velocity through polar coordinates,
    % assuming initial position is at the origin
    dxdt0 = vel0(:, 2).*cos(vel0(:, 1));
    dydt0 = vel0(:, 2).*sin(vel0(:, 1));
    dxdt1 = vel1(:, 2).*cos(vel1(:, 1));
    dydt1 = vel1(:, 2).*sin(vel1(:, 1));
    dxdt2 = vel2(:, 2).*cos(vel2(:, 1));
    dydt2 = vel2(:, 2).*sin(vel2(:, 1));
    dxdt3 = vel3(:, 2).*cos(vel3(:, 1));
    dydt3 = vel3(:, 2).*sin(vel3(:, 1));
    dxdt4 = vel4(:, 2).*cos(vel4(:, 1));
    dydt4 = vel4(:, 2).*sin(vel4(:, 1));
    % finding x and y components of position by integrating
    x0 = cumtrapz(t, dxdt0); % integrate to find x position
    y0 = cumtrapz(t, dydt0); % integrate to find y position
    x1 = cumtrapz(t, dxdt1); % integrate to find x position
    y1 = cumtrapz(t, dydt1); % integrate to find y position
    x2 = cumtrapz(t, dxdt2); % integrate to find x position
    y2 = cumtrapz(t, dydt2); % integrate to find y position
    x3 = cumtrapz(t, dxdt3); % integrate to find x position
    y3 = cumtrapz(t, dydt3); % integrate to find y position
    x4 = cumtrapz(t, dxdt4); % integrate to find x position
    y4 = cumtrapz(t, dydt4); % integrate to find y position
    plot(x0, y0, '-r', 'LineWidth', 2);
    plot(x1, y1, '-g', 'LineWidth', 2);
    plot(x2, y2, '-b', 'LineWidth', 2);
    plot(x3, y3, '-m', 'LineWidth', 2);
    plot(x4, y4, '-k', 'LineWidth', 2);
    xlabel('X Position', 'Interpreter', 'latex');
    ylabel('Y Position', 'Interpreter', 'latex');
    %legend('D = 0','D = 1','D = 2','D = 3','D = 4','Interpreter', 'latex');
    title('Trajectories of Glider for different values of $\theta_0$','Interpreter','latex');
    set(gca, 'FontSize', 44);

figure(1);
xlabel('Time (seconds)', 'Interpreter', 'latex');
ylabel('State Variables', 'Interpreter', 'latex');
legend('show', 'Interpreter', 'latex');
set(gca, 'FontSize', 44);
%set x-axis to constant length so charts are comparable?
title('$\theta$ and $s$ over Time', 'Interpreter', 'latex');

%savegcf to make PDF?