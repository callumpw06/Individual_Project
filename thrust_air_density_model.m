%----- Modelling a Propelled Aeroplane with Altitude-Dependent Thrust -----
close all; clc; clear all;

%----- Defining Parameters -----
D_values = [1 2 3 4];
dt = 0.001;
tspan = 0:dt:15; 
theta_0 = 0;      
H = 3.0;          
k = 1.5;          
T_base = 1.5;     
colors = ['k' 'r' 'g' 'b'];

%----- Setup Figure for Trajectory (Shared) -----
f_traj = figure('Name', 'Trajectory Comparison', 'Color', 'w', 'Theme','light');
ax_traj = gca;
hold(ax_traj, 'on'); grid(ax_traj, 'on');
xlabel(ax_traj, 'Horizontal Position ($x$)', 'Interpreter', 'latex');
ylabel(ax_traj, 'Altitude ($y$)', 'Interpreter', 'latex');
title(ax_traj, ['Trajectory Comparison ($k=' num2str(k) '$, $H=' num2str(H) '$)'], 'Interpreter', 'latex');
set(ax_traj, 'FontSize', 24);
xlim(ax_traj, [0 15]);
ylim(ax_traj, [-1 1.1]);

%----- Loop -----
for i = 1:length(D_values)
    D = D_values(i);
    
    %----- Defining system of Equations -----
    u = @(t,v) [(v(2)^2 * exp(-v(3)/H) - cos(v(1))) / v(2);
        -sin(v(1)) - D*v(2)^2*exp(-v(3)/H) + T_base*exp(-k*v(3)/H);
        v(2) * sin(v(1)) ];
        
    %----- Initial Conditions & Solve -----
    v0 = [theta_0; 1.5; 0]; 
    opts = odeset('RelTol',1e-6,'AbsTol',1e-6);
    [t, vel] = ode45(u, tspan, v0, opts);
    
    figure('Name', ['Dynamics for D = ' num2str(D)], 'Color', 'w', 'Theme', 'light');
    ax_curr = gca; 
    hold(ax_curr, 'on'); grid(ax_curr, 'on');
    set(ax_curr, 'FontSize', 24);
    
    % Labels
    title(ax_curr, ['Flight Variables ($D=' num2str(D) '$)'], 'Interpreter','latex','FontSize', 24);
    xlabel(ax_curr, 'Time, $t$', 'Interpreter', 'latex', 'FontSize', 24);
    ylabel(ax_curr, 'Magnitude', 'Interpreter','latex', 'FontSize', 24);

    % Plotting Variables (All on one axis, using LineStyle to distinguish)
    % 1. Angle (Solid)
    plot(ax_curr, t, vel(:,1), 'LineStyle', '-', 'Color', colors(i), 'LineWidth', 2, ...
        'DisplayName', '$\theta$ (Angle)');
    
    % 2. Speed (Dashed)
    plot(ax_curr, t, vel(:,2), 'LineStyle', '--', 'Color', colors(i), 'LineWidth', 2, ...
        'DisplayName', '$s$ (Speed)');
    
    % 3. Altitude (Dotted)
    plot(ax_curr, t, vel(:,3), 'LineStyle', ':', 'Color', colors(i), 'LineWidth', 2, ...
        'DisplayName', '$y$ (Altitude)');
        
    % 4. Thrust (Dash-Dot)
    plot(ax_curr, t, T_base.*exp(-k.*vel(:,3)/H), "LineStyle", "-." , "Color", colors(i), "LineWidth", 2, ...
        "DisplayName", 'Thrust');
        
    % Add Legend for this specific figure
    legend(ax_curr, 'show', 'Interpreter', 'latex', 'FontSize', 20, 'Location', 'best');

    % Switch focus back to the trajectory figure
    set(0, 'CurrentFigure', f_traj); 
    
    x = cumtrapz(t, vel(:,2).*cos(vel(:,1)));
    plot(ax_traj, x, vel(:,3), 'Color', colors(i), 'LineWidth', 2, ...
        'DisplayName', ['D = ' num2str(D)]);
end

% Final Legend for Trajectory Figure
legend(ax_traj, 'show', 'Interpreter', 'latex', 'FontSize', 24, 'Location', 'southeast');