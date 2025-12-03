clear; clc; close all;

% ----------------------------
% Parameters
% ----------------------------
N     = 5;          % number of agents
dim   = 3;          % 3D positions

kd    = 0.5;        % pursuit / relative position gain
kc    = 0.1;        % stiffness / damping base gain
alpha = pi/N;       % in–plane bearing angle

omega_R = 2*pi/60;  % rad/s, 1 "orbit" = 60 s (conceptual)
kg      = omega_R / (2*sin(pi/N));

R_des = 5;        % desired radius
k_r   = 0.5;      % radial position gain
k_vr  = 0.5;      % radial velocity gain

% Out-of-plane shaping
z_ratio = 0.5;      % y:z scaling (was z0)
phi_z   = pi/4;     % phase between x and z

Tend  = 180;
tspan = [0 Tend];

% ----------------------------
% Initial conditions (3D)
% ----------------------------
x0_pos = 5 * randn(dim, N);
x0_vel = zeros(dim, N);

z0 = [x0_pos(:); x0_vel(:)];   % initial STATE vector

% ----------------------------
% Simulate
% ----------------------------
[t, z] = ode45(@(t,z) cw_cyclic_pursuit_dynamics( ...
                    t, z, N, dim, omega_R, kd, kc, kg, alpha, z_ratio, phi_z), ...
               tspan, z0);

% Extract positions over time
X = z(:, 1:dim*N);

% ----------------------------
% 3D Animation
% ----------------------------
figure; set(gcf,'Color','w');
axis equal; grid on; hold on;
xlabel('x'); ylabel('y'); zlabel('z');
title('CW-style 3D Cyclic Pursuit');

colors = lines(N);
numT   = length(t);

allx = X(:, 1:3:end);
ally = X(:, 2:3:end);
allz = X(:, 3:3:end);
xmin = min(allx(:)); xmax = max(allx(:));
ymin = min(ally(:)); ymax = max(ally(:));
zmin = min(allz(:)); zmax = max(allz(:));
padding = 1;
ax = [xmin-padding xmax+padding ymin-padding ymax+padding zmin-padding zmax+padding];

skip = max(1, floor(numT/400));

for k = 1:skip:numT
    cla; hold on; grid on; axis equal;
    xlabel('x'); ylabel('y'); zlabel('z');
    
    xk = reshape(X(k, :), [dim, N]);
    
    for i = 1:N
        traj = reshape(X(1:k, (i-1)*dim+1:i*dim), [], dim);
        plot3(traj(:,1), traj(:,2), traj(:,3), '--', 'Color', colors(i,:));
        plot3(xk(1,i), xk(2,i), xk(3,i), 'o', ...
              'MarkerFaceColor', colors(i,:), ...
              'MarkerEdgeColor', 'k', ...
              'MarkerSize', 8);
    end
    
    for i = 1:N
        j = mod(i, N) + 1;
        plot3([xk(1,i) xk(1,j)], ...
              [xk(2,i) xk(2,j)], ...
              [xk(3,i) xk(3,j)], '--', ...
              'Color', [0.6 0.6 0.6]);
    end
    
    axis(ax);
    view(35, 25);
    title(sprintf('CW-style 3D Cyclic Pursuit (t = %.2f s)', t(k)));
    drawnow;
end

% ================= Dynamics with CW-style control =================
function dzdt = cw_cyclic_pursuit_dynamics(~, z, N, dim, omega_R, kd, kc, kg, alpha, z_ratio, phi_z)
    r = reshape(z(1:dim*N), [dim, N]);
    v = reshape(z(dim*N+1:end), [dim, N]);
    
    R = [cos(alpha) -sin(alpha) 0;
         sin(alpha)  cos(alpha) 0;
         0           0          1];
    
    T = [1/2                 0            0;
         0                   1            0;
         z_ratio*cos(phi_z)  z_ratio*sin(phi_z) 1];
    
    Tinv = inv(T);
    
    a = zeros(dim, N);
    
    for i = 1:N
        j = mod(i, N) + 1;
        
        ri = r(:,i);   rj = r(:,j);
        vi = v(:,i);   vj = v(:,j);
        
        % CW drift term
        f_i = [ 2*omega_R*vi(2) + 3*omega_R^2*ri(1);
               -2*omega_R*vi(1);
               -omega_R^2*ri(3) ];
        
        % CW-style pursuit control
        pos_term = kd * T*R*Tinv * (rj - ri);
        vel_term =      T*R*Tinv * (vj - vi);
        spring   = -kc*kd        *  ri;
        damp     = -(kc + kd/kg) *  vi;
        
        u_ctrl = -f_i + kg*(pos_term + vel_term + spring + damp);
        
        a(:,i) = f_i + u_ctrl;
    end
    
    drdt = v;
    dvdt = a;
    
    dzdt = [drdt(:); dvdt(:)];
end