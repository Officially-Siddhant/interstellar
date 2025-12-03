clear; clc; close all;

% ----------------------------
% Parameters
% ----------------------------
N     = 5;          % number of agents
dim   = 3;          % 3D positions
kd    = 1.0;        % derivative / pursuit gain
kc    = 0.5;        % "stiffness" gain
alpha = pi/6;       % base rotation angle in radians

Tend  = 40;         % simulation end time
tspan = [0 Tend];

% ----------------------------
% Initial conditions
% ----------------------------
% Random 3D positions (spread out a bit)
x0_pos = 5 * randn(dim, N);

% Initial velocities (can set to zero)
x0_vel = zeros(dim, N);

% Flatten into state vector: [x; v]
z0 = [x0_pos(:); x0_vel(:)];

% ----------------------------
% Simulate
% ----------------------------
[t, z] = ode45(@(t,z) cyclic_pursuit_dynamics(t, z, N, dim, kd, kc, alpha), tspan, z0);

% Extract positions over time
X = z(:, 1:dim*N);              % positions
V = z(:, dim*N+1:end);          % velocities

% ----------------------------
% Animation (3D)
% ----------------------------
figure; hold on;
axis equal;
grid on;
xlabel('x'); ylabel('y'); zlabel('z');
title('3D Cyclic Pursuit with 5 Agents');
set(gcf, 'Color', 'w');

colors = lines(N);
skip   = max(1, floor(length(t)/400));  % frame skip for speed

% Precompute axis limits
allx = X(:, 1:3:end);
ally = X(:, 2:3:end);
allz = X(:, 3:3:end);
xmin = min(allx(:)); xmax = max(allx(:));
ymin = min(ally(:)); ymax = max(ally(:));
zmin = min(allz(:)); zmax = max(allz(:));
padding = 1;
ax = [xmin-padding xmax+padding ymin-padding ymax+padding zmin-padding zmax+padding];

for k = 1:skip:length(t)
    clf; hold on; axis equal; grid on;
    xlabel('x'); ylabel('y'); zlabel('z');
    
    % Current positions (dim x N)
    xk = reshape(X(k, :), [dim, N]);
    
    % Plot agents and traces
    for i = 1:N
        % trajectory up to time t(k)
        traj = reshape(X(1:k, (i-1)*dim+1:i*dim), [], dim);
        plot3(traj(:,1), traj(:,2), traj(:,3), '--', 'Color', colors(i,:)); % trace
        
        % current position
        plot3(xk(1,i), xk(2,i), xk(3,i), 'o', ...
            'MarkerFaceColor', colors(i,:), ...
            'MarkerEdgeColor', 'k', ...
            'MarkerSize', 8);
        text(xk(1,i)+0.1, xk(2,i)+0.1, xk(3,i)+0.1, sprintf('%d', i), ...
             'Color', colors(i,:), 'FontSize', 10);
    end
    
    % Draw lines i -> i+1 (cyclic)
    for i = 1:N
        j = mod(i, N) + 1;
        plot3([xk(1,i) xk(1,j)], ...
              [xk(2,i) xk(2,j)], ...
              [xk(3,i) xk(3,j)], '--', ...
              'Color', [0.6 0.6 0.6]);
    end
    
    axis(ax);
    view(35, 25);
    title(sprintf('3D Cyclic Pursuit (t = %.2f s)', t(k)));
    drawnow;
end

% ================= Dynamics =================
function dzdt = cyclic_pursuit_dynamics(~, z, N, dim, kd, kc, alpha_unused)
    % Unpack state
    x = reshape(z(1:dim*N), [dim, N]);      % positions
    v = reshape(z(dim*N+1:end), [dim, N]);  % velocities
    
    % Parameters
    ka    = 0.5;              % alpha adaptation gain
    R_des = 1.0;              % desired radius (conceptual)
    r_des = 2 * R_des * sin(pi/N);  % desired inter-agent distance (chord)
    
    kv = 2.0;                 % velocity tracking gain

    u = zeros(dim, N);

    for i = 1:N
        j = mod(i, N) + 1;   % agent i+1 (cyclic)

        dx   = x(:,j) - x(:,i);
        d_ij = norm(dx);

        % Adaptive in-plane bearing angle (scalar)
        alpha_i = (pi * (1/N)) + ka * (r_des - d_ij);

        % 3D rotation about z-axis
        R = [cos(alpha_i) -sin(alpha_i) 0;
             sin(alpha_i)  cos(alpha_i) 0;
             0             0            1];

        % --------- SINGLE-INTEGRATOR FIELD (desired velocity) ----------
        v_des = kd * R * dx - kc * x(:,i);

        % --------- DOUBLE-INTEGRATOR EXTENSION ------------------------
        % Make acceleration pull velocity toward v_des
        u(:,i) = -kv * (v(:,i) - v_des);
    end

    % Double integrator dynamics
    xdot = v;
    vdot = u;

    dzdt = [xdot(:); vdot(:)];
end