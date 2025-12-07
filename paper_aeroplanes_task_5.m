close all;
%set D for plot
%PLOTS FOR D 1 TO 4, EQUILIBRIUM POINTS
D = [0, 1, 2, 3, 4];
%creates coordinate grid
%adjust grid size?
[X, Y, Z] = meshgrid(-pi:pi/20:pi,0:0.05:2,-pi:pi/20:pi);
%remove values where X (theta) does not equal Z in order to create the helix
%pattern for the phase cylinder
Xtemp = X;
Ztemp = Z;
X(Xtemp~=Ztemp) = 0;
Z(Xtemp~=Ztemp) = 0;
%convert to cylindrical coordinates
Xpolr = Y.*cos(X);
Ypolr = Y.*sin(X);
%separate grid for calculating trajectory
%adjust grid size?
[Xstream, Ystream, Zstream] = meshgrid(-pi:pi/20:pi,0:0.05:2,-pi:pi/20:pi);
%vectors setup for trajectory calculation
[theta_stream, s_stream] = cart2pol(Xstream,Ystream);
%equilibrium points for each value of D, calculated analytically beforehand
equilibrium_theta = [atan(0),atan(-1),atan(-2),atan(-3),atan(-4)];
equilibrium_s = sqrt(cos(equilibrium_theta));
equilibrium_x = equilibrium_s.*cos(equilibrium_theta);
equilibrium_y = equilibrium_s.*sin(equilibrium_theta);
%plot breaks with scaling when step divides domain
%scale factor for plotting DOESNT DO ANYTHING, FIX TO MAKE LESS CLUTTERED
sf = 0.1;
for i = 1:length(D)
    figure('Theme', 'light');
    %set up functions from X,Y,Z -> U,V,W
    U = (Y.^2 - cos(X))./Y; %normalise vectors to length 1, scaling afterwards. show length with colours?
    V = -sin(X) - D(i)*Y.^2;
    %don't want the vectors to move vertically
    W = zeros(size(U));
    Upolr = V.*cos(U);
    Vpolr = V.*sin(U);
    %setup grids for U direction, V direction, and Magnitude
    magnitude = hypot(Upolr,Vpolr);
    %Scale U and V
    Unorm = sf.*Upolr./magnitude;
    Vnorm = sf.*Vpolr./magnitude;
    hold on;
    %trajectory calculation
    [theta_stream, s_stream] = cart2pol(Xstream,Ystream);
    Ustream = (s_stream.^2 - cos(theta_stream))./s_stream;
    Vstream = -sin(theta_stream) - D(i)*s_stream.^2;
    %how to make W follow the helix?
    Wstream = (s_stream.^2 - cos(Zstream))./s_stream;
    %Plot
    quiver3(Xpolr,Ypolr,Z,Unorm,Vnorm,W,0, 'HandleVisibility', 'off', 'Color', 'black');
    %shows how system develops from initial conditions
    %ADJUST INITIAL CONDITIONS TO MATCH CYLINDRICAL COORDS INSTEAD OF
    %CARTESIAN COORDS
    streamline(Xstream,Ystream,Zstream,Ustream,Vstream,Wstream,cos(0),sin(0),0, 'LineWidth', 2, 'DisplayName', '$(s_0,\theta_0)=(1,0)$', 'Color', 'r');
    streamline(Xstream,Ystream,Zstream,Ustream,Vstream,Wstream,cos(pi/8),sin(pi/8),pi/8, 'LineWidth', 2, 'DisplayName', '$(s_0,\theta_0)=(1,\pi/8)$', 'Color', 'g');
    streamline(Xstream,Ystream,Zstream,Ustream,Vstream,Wstream,cos(pi/4),sin(pi/4),pi/4, 'LineWidth', 2, 'DisplayName', '$(s_0,\theta_0)=(1,\pi/4)$', 'Color', 'b');
    %plot contour line and/or equilibrium points with scatter
    plot3(equilibrium_x(i),equilibrium_y(i),equilibrium_theta(i), '.', 'Color', 'm', 'LineWidth', 2, 'HandleVisibility', 'off','DisplayName','Equilibrium Point');
    %RELABEL AXES, POSSIBLY MAKE AXES A CYLINDER SHAPE?
    %axis equal;
    %xlabel('$\theta$', 'Interpreter', 'latex');
    %ylabel('$s$', 'Interpreter', 'latex', 'Rotation', 0);
    zlabel('$z=\theta$', 'Interpreter', 'latex');
    set(gca, 'FontSize', 42);
    %EXPAND AXES? MIGHT BE CUTTING OFF THE HEADS OF SOME VECTORS
    axis([-2 2 -2 2 -pi pi]);
    legend('show','Location', 'best', 'Interpreter', 'latex');
    title(sprintf('Phase Plane ($D = %i$)', i-1), 'Interpreter', 'latex');
end