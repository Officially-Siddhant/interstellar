clear; clc; close all;

% ----------------------------
% Parameters
% ----------------------------
N     = 5;          % number of agents
dim   = 2;          % 2D positions
kd    = 1.0;        % derivative / pursuit gain
kc    = 0.5;        % "CW-style" stiffness gain
alpha = pi/6;       % rotation angle in radians

Tend  = 40;         % simulation end time
tspan = [0 Tend];

% ----------------------------
% Initial conditions
% ----------------------------
% Random 2D positions (spread out a bit)
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
% Animation
% ----------------------------
figure; hold on;
axis equal;
grid on;
xlabel('x'); ylabel('y');
title('Cyclic Pursuit with 5 Agents');
set(gcf, 'Color', 'w');

colors = lines(N);
skip   = max(1, floor(length(t)/400));  % frame skip for speed

for k = 1:skip:length(t)
    clf; hold on; axis equal;
    
    % Current positions (dim x N)
    xk = reshape(X(k, :), [dim, N]);
    
    % Plot agents
    for i = 1:N
        plot(xk(1,i), xk(2,i), 'o', ...
            'MarkerFaceColor', colors(i,:), ...
            'MarkerEdgeColor', 'k', ...
            'MarkerSize', 8);
        text(xk(1,i)+0.1, xk(2,i)+0.1, sprintf('%d', i), ...
             'Color', colors(i,:), 'FontSize', 10);
    end
    
    % Draw lines i -> i+1 (cyclic)
    for i = 1:N
        j = mod(i, N) + 1;
        plot([xk(1,i) xk(1,j)], [xk(2,i) xk(2,j)], '--', ...
             'Color', [0.6 0.6 0.6]);
    end
    
    % Set axes
    allx = X(:, 1:2:end);  % x-coordinates of all agents across time
    ally = X(:, 2:2:end);  % y-coordinates of all agents across time
    xmin = min(allx(:)); xmax = max(allx(:));
    ymin = min(ally(:)); ymax = max(ally(:));
    padding = 1;
    axis([xmin-padding xmax+padding ymin-padding ymax+padding]);
    
    title(sprintf('Cyclic Pursuit (t = %.2f s)', t(k)));
    drawnow;
end

function dzdt = cyclic_pursuit_dynamics(~, z, N, dim, kd, kc, alpha_unused)
    % Unpack state
    x = reshape(z(1:dim*N), [dim, N]);      % positions
    v = reshape(z(dim*N+1:end), [dim, N]);  % velocities
    
    % Parameters
    ka = 0.5;              % your alpha adaptation gain
    R_des = 1.0;           % desired circle radius (if you want)
    r = 2 * R_des * sin(pi/N);  % desired inter-agent distance (chord)
    
    kv = 2.0;              % velocity tracking gain (tune this)

    u = zeros(dim, N);

    for i = 1:N
        j = mod(i, N) + 1;   % agent i+1 (cyclic)

        dx = x(:,j) - x(:,i);
        d_ij = norm(dx);

        % Your ever-changing alpha, but scalar and dimensionally OK
        alpha_i = (pi * (1/N)) + ka * (r - d_ij);

        % Rotation matrix with that alpha
        R = [cos(alpha_i) -sin(alpha_i);
             sin(alpha_i)  cos(alpha_i)];

        % --------- SINGLE-INTEGRATOR FIELD (desired velocity) ----------
        % You can pack all your "CW-style" stuff here if you like
        v_des = kd * R * dx - kc * x(:,i);   % attraction toward the origin / radius control

        % --------- DOUBLE-INTEGRATOR EXTENSION ------------------------
        % Make acceleration pull velocity toward v_des
        u(:,i) = -kv * (v(:,i) - v_des);
    end

    % Double integrator dynamics
    xdot = v;
    vdot = u;

    dzdt = [xdot(:); vdot(:)];
end