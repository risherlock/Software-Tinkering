%{
  Reproduction of Peter Lynch's implementation.
  peniac: Peter's eniac

  References:
    [1] Numerical Integration of the Barotropic Vorticity Equation
        by Charney, Fjortoft, Neumann (1950)
    [2] Peter Lynch - The ENIAC Forecasts: A Recreation (2008)

  2023-12-20
%}

% clc
clear
close all

% Simulation time params
t = 86400; % Simulation time, s
dt = 3600; % Time step, s

% Simulation grid params
nx = 19;
ny = 16;
px = 10;
py = 14;
dp = 736e3;
ctrl = -90;

% Initial height
z0 = load('data/1_1949010503.z00');

% Physical params
g = 9.80665; % Acc. due to gravity, m / s^2
a = (4 * 10^7) / (2 * pi); % Radius of Earth, km
w = 2 * pi / (24 * 60 * 60); % Earth's angular rate, r/s

% Variables to be computed for each grid
lat = zeros(nx, ny);
lon = zeros(nx, ny);
m = zeros(nx, ny);
f = zeros(nx, ny);
h = zeros(nx, ny);
eig = zeros(nx-2, ny-2);
sx = zeros(nx-2, nx-2);
sy = zeros(ny-2, ny-2);

% Computations for each grid
for j = 1 : ny
  for i = 1: nx
    % Latitude and longitude
    x = (i - px) * dp;
    y = (j - py) * dp;
    r = sqrt(x^2 + y^2);
    theta = atan2(y, x);
    lambda = theta + deg2rad(90 + ctrl);

    if (lambda > pi)
      lambda = lambda - 2 * pi;
    end

    phi = 2 * ((pi/4) - atan(r / (2 * a)));
    lat(i, j) = lambda;
    lon(i, j) = phi;

    % Eqn. (15)
    m(i, j) = 2 / (1 + sin(phi));
    f(i, j) = 2 * w * sin(phi);
    h(i, j) = g * m(i, j)^2 / f(i, j);

    % Eigenvalues of Laplacian operator
    if ((i < nx-1) && (j < ny-1))
      eig(i, j) = (-4 / dp^2) ...
                * ((sin(0.5 * pi * i / (nx - 1)))^2 ...
                +  (sin(0.5 * pi * j / (ny - 1)))^2);
    end
  end
end

% Parameters for Poission solver: sx, sy
sx = zeros(nx-2, nx-2);
for i = 1 : nx-2
  for j = 1 : nx-2
    sx(i, j) = sin(pi * i * j / (nx - 1));
  end
end

sy = zeros(ny-2, ny-2);
for i = 1 : ny-2
  for j = 1 : ny-2
    sy(i, j) = sin(pi * i * j / (ny - 1));
  end
end

% Variables
xi = zeros(nx, ny);
xi_dot = zeros(nx, ny);
xi_last = zeros(nx, ny);

z = zeros(nx, ny);
z_dot = zeros(nx, ny);
z_xdot = zeros(nx, ny);
z_ydot = zeros(nx, ny);
z_last = zeros(nx, ny);

eta  = zeros(nx, ny);
eta_dot = zeros(nx, ny);
eta_xdot = zeros(nx, ny);
eta_ydot = zeros(nx, ny);

z = z0;
dt_temp = dt;
time = 0 : dt : t;

% Numerical integration loop
for i = 1 : length(time)
  % First and second derivatives of z wrt. x and y
  z_xdot(2:nx-1, :) = (z(3:nx, :) - z(1:nx-2, :)) / (2 * dp);
  z_ydot(:, 2:ny-1) = (z(:, 3:ny) - z(:, 1:ny-2)) / (2 * dp);
  z_xddot(2:nx-1, :) = (z(3:nx, :) + z(1:nx-2, :) - 2 * z(2:nx-1, :)) / dp^2;
  z_yddot(:, 2:ny-1) = (z(:, 3:ny) + z(:, 1:ny-2) - 2 * z(:, 2:ny-1)) / dp^2;

  % Laplacian of height, xi
  xi(2:nx-1, 2:ny-1) = z_xddot(2:nx-1, 2:ny-1) + z_yddot(2:nx-1, 2:ny-1);

  % Set initial condition
  if (i == 1)
    dt = dt_temp / 2;

    % Boundary condition
    xi0 = xi; % Extract Laplacians
    xi0(nx, :) = 2 * xi0(nx-1, :) - xi0(nx-2, :); % E
    xi0(1, :) = 2 * xi0(2, :) - xi0(3, :); % W
    xi0(:, ny) = 2 * xi0(:, ny-1) - xi0(:, ny-2); %  N
    xi0(:, 1) = 2 * xi0(:, 2) - xi0(:, 3); % S
    xi = xi0; % Set BCs

    xi_last = xi0;
    z_last = z0;
    continue;
  end

  % eta and its derivatives
  eta = h .* xi + f;
  eta_xdot(2:nx-1, :) = (eta(3:nx, :) - eta(1:nx-2, :)) / (2 * dp);
  eta_ydot(:, 2:ny-1) = (eta(:, 3:ny) - eta(:, 1:ny-2)) / (2 * dp);

  % xi_dot
  J = eta_xdot .* z_ydot - eta_ydot .* z_xdot;
  xi_dot_tmp = J(2:nx-1, 2:ny-1);

  % Solve Poisson equation using Fourier transform
  fxi_dot = sx * xi_dot_tmp * sy; % FT of xi_dot
  fz_dot = fxi_dot ./ eig; % FT of z_dot

  % Inv FT of fz_dot
  z_dot_tmp = (4 / ((nx - 1) * (ny - 1))) * sx * fz_dot * sy;

  z_dot(2:nx-1, 2:ny-1) = z_dot_tmp;
  xi_dot(2:nx-1, 2:ny-1) = xi_dot_tmp;

  % Compute xi_dot at boundaries
  sign_e = max(0, -sign(z_ydot(nx,:)));
  sign_w = max(0, +sign(z_ydot(1,:)));
  sign_n = max(0, +sign(z_xdot(:,ny)));
  sign_s = max(0, -sign(z_xdot(:,1)));
  xi_dot(nx,:) = sign_e .* (2 * xi_dot(nx-1, :) - xi_dot(nx-2, :)); % E
  xi_dot(1,:)  = sign_w .* (2 * xi_dot(2,:) - xi_dot(3,:)); % W
  xi_dot(:,ny) = sign_n .* (2 * xi_dot(:, ny-1) - xi_dot(:, ny-2)); % N
  xi_dot(:,1)  = sign_s .* (2 * xi_dot(:,2) - xi_dot(:, 3)); % S

  % Inegration
  z_ans = z_last +  2 * dt * z_dot;
  xi_ans = xi_last + 2 * dt * xi_dot;

  % Memory
  z_last = z;
  xi_last = xi;
  z = z_ans;
  xi = xi_ans;

  % Restore dt
  dt = dt_temp;
end

fprintf("Forcast complete!\n");
