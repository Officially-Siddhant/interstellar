%% CW-style 3D cyclic pursuit with radius-keeping
clear; clc; close all;

% ----------------------------
% Parameters
% ----------------------------
N     = 5;          % number of agents
dim   = 3;          % 3D positions

kd    = 0.5;        % pursuit / relative position gain
kc    = 0.1;        % base damping gain
alpha = pi/N;       % in–plane bearing angle

% "Orbital" rate (scaled: 1 orbit ~ 60 s)
omega_R = 2*pi/60;          % [rad/s]
kg      = omega_R / (2*sin(pi/N));

% Out-of-plane shaping
z_ratio = 0.5;              % y:z scaling
phi_z   = pi/4;             % phase between x and z

% Desired radius and radial PD gains
R_des = 5.0;                % desired radius from origin
k_r   = 0.5;                % radial position gain
k_vr  = 0.5;                % radial velocity gain

% Time
Tend  = 60;                 % simulate 60 s
tspan = [0 Tend];

% ----------------------------
% Initial conditions (3D)
% ----------------------------
r0 = 5 * randn(dim, N);     % random 3D positions
v0 = zeros(dim, N);         % initial velocities

z0 = [r0(:); v0(:)];        % state vector [r; v]

% ----------------------------
% Simulate
% ----------------------------
[t, z] = ode45(@(t,z) cw_cyclic_pursuit_radius( ...
                    t, z, N, dim, omega_R, kd, kc, kg, alpha, ...
                    z_ratio, phi_z, R_des, k_r, k_vr), ...
               tspan, z0);

% Extract positions over time
X = z(:, 1:dim*N);          % positions (flattened)

% ----------------------------
% 3D Animation
% ----------------------------
figure; set(gcf,'Color','w');
axis equal; grid on; hold on;
xlabel('x'); ylabel('y'); zlabel('z');
title('CW-style 3D Cyclic Pursuit with Radius Keeping');

colors = lines(N);
numT   = length(t);

% Precompute axis limits
allx = X(:, 1:3:end);
ally = X(:, 2:3:end);
allz = X(:, 3:3:end);
xmin = min(allx(:)); xmax = max(allx(:));
ymin = min(ally(:)); ymax = max(ally(:));
zmin = min(allz(:)); zmax = max(allz(:));
padding = 1;
ax = [xmin-padding xmax+padding ...
      ymin-padding ymax+padding ...
      zmin-padding zmax+padding];

skip = max(1, floor(numT/400));
playback_scale = 1.0;   % 1.0 ≈ real-time, >1 faster, <1 slower

for k = 1:skip:numT
    cla; hold on; grid on; axis equal;
    xlabel('x'); ylabel('y'); zlabel('z');
    
    rk = reshape(X(k, :), [dim, N]);
    
    for i = 1:N
        traj = reshape(X(1:k, (i-1)*dim+1:i*dim), [], dim);
        plot3(traj(:,1), traj(:,2), traj(:,3), '--', 'Color', colors(i,:)); % trace
        
        plot3(rk(1,i), rk(2,i), rk(3,i), 'o', ...
              'MarkerFaceColor', colors(i,:), ...
              'MarkerEdgeColor', 'k', ...
              'MarkerSize', 8);
    end
    
    for i = 1:N
        j = mod(i, N) + 1;
        plot3([rk(1,i) rk(1,j)], ...
              [rk(2,i) rk(2,j)], ...
              [rk(3,i) rk(3,j)], '--', ...
              'Color', [0.6 0.6 0.6]);
    end
    
    axis(ax);
    view(35, 25);
    title(sprintf('CW-style 3D Cyclic Pursuit (t = %.2f s)', t(k)));
    drawnow;
    
    if k > 1
        dt_sim = t(k) - t(k-1);
        pause(dt_sim / playback_scale);
    end
end


%% ================= Dynamics with CW-style pursuit + radius PD =================
function dzdt = cw_cyclic_pursuit_radius(~, z, N, dim, omega_R, kd, kc, kg, alpha, ...
                                         z_ratio, phi_z, R_des, k_r, k_vr)
    % Unpack state
    r = reshape(z(1:dim*N), [dim, N]);      % positions (x,y,z)
    v = reshape(z(dim*N+1:end), [dim, N]);  % velocities
    
    % Rotation and shaping matrices
    R = [cos(alpha) -sin(alpha) 0;
         sin(alpha)  cos(alpha) 0;
         0           0          1];
    
    T = [1/2                 0                  0;
         0                   1                  0;
         z_ratio*cos(phi_z)  z_ratio*sin(phi_z) 1];
    
    Tinv = inv(T);
    
    a = zeros(dim, N);   % accelerations
    
    for i = 1:N
        j = mod(i, N) + 1;   % cyclic neighbor
        
        ri = r(:,i);   rj = r(:,j);
        vi = v(:,i);   vj = v(:,j);
        
        % --- CW drift term f_i ---
        f_i = [ 2*omega_R*vi(2) + 3*omega_R^2*ri(1);
               -2*omega_R*vi(1);
               -omega_R^2*ri(3) ];
        
        % --- Cyclic pursuit (relative position + velocity) ---
        pos_term = kd * T*R*Tinv * (rj - ri);
        vel_term =      T*R*Tinv * (vj - vi);
        
        % --- Radial PD to keep ||ri|| ≈ R_des ---
        r_norm = norm(ri) + 1e-9;        % avoid divide by zero
        n_i    = ri / r_norm;            % radial unit vector
        e_r    = r_norm - R_des;         % radius error
        v_r    = dot(vi, n_i);           % radial velocity
        
        radial_PD = -k_r*e_r*n_i - k_vr*v_r*n_i;
        
        % --- Global damping (optional) ---
        damp = -(kc + kd/kg) * vi;
        
        % Control law: u_ctrl acts in CW frame
        u_ctrl = -f_i + kg*(pos_term + vel_term + radial_PD + damp);
        
        % Total acceleration r̈ = f_i + u_ctrl
        a(:,i) = f_i + u_ctrl;
    end
    
    drdt = v;
    dvdt = a;
    
    dzdt = [drdt(:); dvdt(:)];
end