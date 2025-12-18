%----- Modelling Supersonic Flight: Mach Tuck & Ackeret's Theory -----
close all; clc; clear all;

%==========================================================================
%   PARAMETERS & CONSTANTS
%==========================================================================
% Sources: 
% [1] Anderson, "Fundamentals of Aerodynamics" (Speed of Sound, Ackeret)
% [2] Raymer, "Aircraft Design" (Drag Divergence, CP Shift)

% Environmental
s_v = 340;          % Speed of Sound (m/s) [Standard Atmosphere]

% Aerodynamic Coefficients
D_sub = 0.02;       % Zero-lift drag (Subsonic - roughly F-16 clean)
D_super = 0.05;     % Wave drag penalty (Supersonic)
k_drag = 0.5;       % Width of the transonic drag rise

% Aircraft Specs
m = 10000;          % Mass (kg)
S_wing = 30;        % Wing Area (m^2)
g = 9.81;           % Gravity (m/s^2)

% Engine Specs (Afterburner)
T_max = 160000;     % Max Thrust (N) - High T/W ratio to break barrier
k_thrust = 1.0e-4;  % Thrust lapse with altitude

% Simulation Control
dt = 0.01;
tspan = 0:dt:20;    

% PILOT CONTROL (Trim)
% Pilot holds nose 10 degrees above horizon
pitch_attitude = 10 * (pi/180); 

%==========================================================================
%   SOLVER
%==========================================================================
% Initial Conditions: [Theta=0, Velocity=250m/s, Altitude=1000m]
v0 = [0; 250; 1000]; 

% Pass parameters to dynamics function
dynamics = @(t, v) plane_dynamics(v, s_v, m, S_wing, g, T_max, k_thrust, ...
                                  D_sub, D_super, k_drag, pitch_attitude);

opts = odeset('RelTol', 1e-6, 'AbsTol', 1e-6);
[t, state] = ode45(dynamics, tspan, v0, opts);

%==========================================================================
%   DATA EXTRACTION & PLOTTING
%==========================================================================
theta = state(:,1);
v = state(:,2);
h = state(:,3);
Mach = v ./ s_v;

% Re-calculate Alpha and Lift for visualization
alpha = pitch_attitude - theta; 
L_vals = zeros(size(t));
for i = 1:length(t)
    [~, L_vals(i)] = get_aero_forces(theta(i), v(i), h(i), ...
                     s_v, m, S_wing, D_sub, D_super, k_drag, pitch_attitude);
end

% --- PLOT SETUP ---
figure('Name', 'Supersonic Dynamics Model', 'Color', 'w', 'Position', [100, 100, 1200, 800], 'Theme', 'light');

% 1. Velocity & Mach
subplot(2,2,1); hold on; grid on;
plot(t, Mach, 'b', 'LineWidth', 2);
yline(1.0, 'r--', 'Sound Barrier (Mach 1)');
xlabel('Time (s)'); ylabel('Mach Number');
title('Acceleration Phase');
set(gca, 'FontSize', 12);

% 2. Flight Path (The Mach Tuck)
subplot(2,2,2); hold on; grid on;
plot(t, theta*(180/pi), 'k', 'LineWidth', 2);
yline(pitch_attitude*(180/pi), 'g:', 'Pilot Stick Input (Pitch)');
xline(t(find(Mach>1,1)), 'r--', 'Transition');
xlabel('Time (s)'); ylabel('Flight Path Angle \theta (deg)');
title('Flight Path Stability (Note the Dip!)');
legend('Path', 'Pitch Input', 'Mach 1', 'Location', 'SouthEast');
set(gca, 'FontSize', 12);

% 3. Lift Force
subplot(2,2,3); hold on; grid on;
plot(t, L_vals, 'b', 'LineWidth', 2);
yline(m*g, 'k--', 'Aircraft Weight');
xlabel('Time (s)'); ylabel('Lift Force (N)');
title('Lift vs Weight');
set(gca, 'FontSize', 12);

% 4. Altitude
subplot(2,2,4); hold on; grid on;
plot(t, h, 'm', 'LineWidth', 2);
xlabel('Time (s)'); ylabel('Altitude (m)');
title('Altitude Profile');
set(gca, 'FontSize', 12);


%==========================================================================
%   DYNAMICS FUNCTIONS
%==========================================================================
function dv = plane_dynamics(v_state, s_v, m, S, g, T_max, k_thrust, ...
                             D_sub, D_super, k_drag, pitch_attitude)
    theta = v_state(1);
    v = v_state(2);
    h = v_state(3);
    
    % Get Aerodynamic Forces (Lift & Drag)
    [Drag, Lift] = get_aero_forces(theta, v, h, s_v, m, S, ...
                                   D_sub, D_super, k_drag, pitch_attitude);
    
    % Engine Thrust (Decreases with altitude)
    Thrust = T_max * exp(-k_thrust * h);
    
    % --- EQUATIONS OF MOTION ---
    % 1. Flight Path Angle (dTheta/dt)
    % Force Balance: (Lift + Thrust*sin(alpha) - Weight*cos(theta)) / mv
    % Note: We assume small alpha, so Thrust acts roughly along path
    d_theta = (Lift - m*g*cos(theta)) / (m*v);
    
    % 2. Acceleration (dv/dt)
    % Force Balance: (Thrust - Drag - Weight*sin(theta)) / m
    d_v = (Thrust - Drag - m*g*sin(theta)) / m;
    
    % 3. Climb Rate (dh/dt)
    d_h = v * sin(theta);
    
    dv = [d_theta; d_v; d_h];
end

function [Drag, Lift] = get_aero_forces(theta, v, h, s_v, m, S, ...
                                        D_sub, D_super, k_drag, pitch_attitude)
    rho = 1.225 * exp(-h/8500); % Air Density (Exponential Atmos)
    q = 0.5 * rho * v^2;        % Dynamic Pressure
    M = v / s_v;                % Mach Number
    
    % Calculate Angle of Attack (Alpha)
    alpha = pitch_attitude - theta; 
    
    % --- DRAG MODEL (Sigmoid Rise) ---
    % CD rises from 0.02 to 0.05 as we cross Mach 1
    Cd = D_sub + (D_super - D_sub) / (1 + exp(-k_drag * (v - s_v)));
    % Add induced drag (drag due to lift)
    k_induced = 0.1;
    Cd_total = Cd + k_induced * alpha^2;
    Drag = Cd_total * q * S;
    
    % --- LIFT MODEL (The Core Physics) ---
    if M < 1.0
        % [SUBSONIC]
        % Linear Lift Slope: CL = 2*pi * alpha
        Cl = 2 * pi * alpha;
        
    else
        % [SUPERSONIC] - Ackeret's Theory
        % CL = 4 * alpha / sqrt(M^2 - 1)
        
        term = sqrt(M^2 - 1);
        % Mathematical safety clamp (Prandtl-Glauert Singularity)
        if term < 0.3, term = 0.3; end 
        
        Cl = (4 * alpha) / term;
        
        % [MACH TUCK EFFECT]
        % In a real plane, the Center of Pressure moves back. 
        % This reduces the pitching moment. In a point-mass model, 
        % we simulate this as a reduction in effective Lift Coefficient
        % or a penalty requiring higher alpha.
        
        Tuck_Severity = 0.5; % Tunable parameter for "Stability"
        Cl = Cl - Tuck_Severity * (M - 1.0); 
    end
    
    Lift = Cl * q * S;
end