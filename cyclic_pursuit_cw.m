% Clohessy Wilshire Model

% Steps to take
% 1. define system parameters such as orbital duration, gain vectors for
% the control law being used, N = 5 agents

% Define initial start conditions X_0 = zeros(3,N)
% Define simulation duration 
% μ = gravitational constant
% Define omega_R = andular velocity
% Define radius R_ref = (μ / omega^2)^(1/3)
% simulate the system by calling a function like [t, z] = ode45(@(t,z)
% cyclic_pursuit_dynamics(t, z, N, dim, kd, kc, alpha), tspan, z0); that
% solves the dynamics based on the model provide

% control law is given and will be defined inside the function

% animate - earth = 3D sphere

% function dXdt = cw_dynamics(t, X, N, dim, kd, kc, alpha_unused)
%     x = X(1,N); % x components
%     y = X(2,N); % y components
%     z = X(3,N); % z components
%     dXdt(1) = ...;
%     dXdt(2) = ...;
%     dXdt(3) = ...;
% end

%% CW cyclic pursuit with 60 s scaled orbit
clear; clc; close all;

% ----------------------------
% System / formation params
% ----------------------------
N = 5;
dim = 3;          % 3D

% Scale: 1 orbit = 60 s means that our omega_R = 2*pi / 60
omega_R = 2*pi/60;

% Controller gains
kd    = 0.5;
kc    = 0.1;
alpha = pi/N;

% Out-of-plane shaping parameters
z0    = 0.5; 
phi_z = pi/4;

% Cyclic pursuit gain constant
kg = omega_R / (2*sin(pi/N));

% ----------------------------
% Initial conditions
% ----------------------------
X0_pos = 100*randn(dim, N);   % relative positions [arbitrary units]
X0_vel = zeros(dim, N);       % initial velocities
X0 = [X0_pos(:); X0_vel(:)];  % 6N×1 state vector

% ----------------------------
% Simulate one "orbit" (60 s)
% ----------------------------
Tend  = 180;
tspan = [0 Tend];

[t, X] = ode45(@(t,X) cw_dynamics(t, X, N, omega_R, kd, kc, kg, alpha, z0, phi_z), ...
               tspan, X0);

% Reshaping positions for plotting
numT = length(t);
% [time, xyz, agent]
r_all = reshape(X(:, 1:dim*N), [numT, dim, N]);   


% ----------
% Animation
% ----------
figure; set(gcf,'Color','w');
axis equal; grid on; hold on;
xlabel('x'); ylabel('y'); zlabel('z');
title('CW cyclic pursuit (scaled, 1 orbit = 60 s)');
colors = lines(N);

% Axis limits
allx = r_all(:,1,:);
ally = r_all(:,2,:);
allz = r_all(:,3,:);
xmin = min(allx(:)); xmax = max(allx(:));
ymin = min(ally(:)); ymax = max(ally(:));
zmin = min(allz(:)); zmax = max(allz(:));
pad = 0.1 * max([xmax-xmin, ymax-ymin, zmax-zmin]);

numT = length(t);
skip = 1;
Tend = t(end);
playback_scale = 1.5;

for k = 1:skip:numT
    cla; hold on;
    for i = 1:N
        traj = squeeze(r_all(1:k, :, i));  % [k x 3]
        plot3(traj(:,1), traj(:,2), traj(:,3), '--', ...
              'Color', colors(i,:));       % dashed trace
        plot3(traj(end,1), traj(end,2), traj(end,3), 'o', ...
              'MarkerFaceColor', colors(i,:), ...
              'MarkerEdgeColor', 'k', ...
              'MarkerSize', 6);
    end
    
    axis([xmin-pad xmax+pad ymin-pad ymax+pad zmin-pad zmax+pad]);
    view(35, 25);
    title(sprintf('CW cyclic pursuit, t = %.1f s', t(k)));
    drawnow;

    % make animation last about Tend / playback_scale seconds
    if k > 1
        dt_sim  = t(k) - t(k-1);          % simulated time step
        pause(dt_sim / playback_scale);   % wall-clock pause
    end
end


%% ===================== Dynamics + control ============================
function dXdt = cw_dynamics(~, X, N, omega_R, kd, kc, kg, alpha, z0, phi_z)
    dim = 3;

    % Unpack positions and velocities: X = [r; v]
    r = reshape(X(1:dim*N), [dim, N]);      % positions (x,y,z) of all agents
    v = reshape(X(dim*N+1:end), [dim, N]);  % velocities

    % Rotation + shaping matrices
    R = [cos(alpha) -sin(alpha) 0;
         sin(alpha)  cos(alpha) 0;
         0           0          1];

    T = [1/2           0          0;
         0             1          0;
         z0*cos(phi_z) z0*sin(phi_z) 1];

    Tinv = inv(T);

    dvdt = zeros(dim, N);

    for i = 1:N
        j = mod(i, N) + 1;  % i+1 (cyclic)

        ri = r(:,i);   rj = r(:,j);
        vi = v(:,i);   vj = v(:,j);

        % -------- HCW drift term f(r_i, v_i) --------
        f_i = [ 2*omega_R*vi(2) + 3*omega_R^2*ri(1);   % x¨
               -2*omega_R*vi(1);                       % y¨
               -omega_R^2*ri(3) ];                     % z¨

        % -------- Cyclic pursuit control law (u_ctrl) --------
        pos_term = kd * T*R*Tinv*(rj - ri);
        vel_term =      T*R*Tinv*(vj - vi);
        spring   = -kc*kd * ri;
        damp     = -(kc + kd/kg) * vi;

        u_ctrl = -f_i + kg*(pos_term + vel_term + spring + damp);

        % Total acceleration: r̈ = f + u_ctrl
        dvdt(:,i) = f_i + u_ctrl;
    end

    drdt = v;
    dXdt = [drdt(:); dvdt(:)];
end