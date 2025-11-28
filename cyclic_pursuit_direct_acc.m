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
V = z(:, dim*N+1:end);          % velocities (unused for plotting here)

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

%% Dynamics function
function dzdt = cyclic_pursuit_dynamics(~, z, N, dim, kd, kc, alpha)
    % Unpack state
    x = reshape(z(1:dim*N), [dim, N]);      % positions
    v = reshape(z(dim*N+1:end), [dim, N]);  % velocities
    
    % Parameters for radius / spacing
    ka = 1.0;              % your gain
    R_des = 1.0;           % desired circle radius
    r = 2 * sin(pi/N);     % inter-agent distance
    
    % r = 2 * R_des * sin(pi/N);
    
    R = [cos(alpha) -sin(alpha);
             sin(alpha)  cos(alpha)];

    % Control for each agent
    u = zeros(dim, N);
    for i = 1:N
        j = mod(i, N) + 1;  % agent i+1 (cyclic)
        
        dx = x(:,j) - x(:,i);
        dv = v(:,j) - v(:,i);
        
        % Actual inter-agent distance (scalar)
        d_ij = norm(dx);
        alpha_i = (pi * (1/N)) + ka * (r - d_ij);
        
        % Rotation matrix
        R2 = [cos(alpha_i) -sin(alpha_i);
             sin(alpha_i)  cos(alpha_i)];
        
        % Given control law:
        % u_i = k_d * R(α_i)(x_{i+1} - x_i)
        %     +      R(α_i)(v_{i+1} - v_i)
        %     - (k_c * k_d) * x_i
        %     - (k_c + k_d) * v_i
        u(:,i) = kd * R * dx ...
               +      R2 * dv ...
               - (kc*kd)   * x(:,i) ...
               - (kc + kd) * v(:,i);
    end
    
    % Double integrator dynamics:
    xdot = v;
    vdot = u;
    
    dzdt = [xdot(:); vdot(:)];
end