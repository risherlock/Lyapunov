% Numerical simulation of quadrotor dynamics with PD control
% 2023-12-02

clc
clear
close all

function w = pd_control(phy_params, state_d, state, K)
  % Physical params
  g = phy_params(1);
  m = phy_params(2);
  l = phy_params(3);
  k = phy_params(4);
  b = phy_params(5);
  Im = phy_params(6);
  Ixx = phy_params(7);
  Iyy = phy_params(8);
  Izz = phy_params(9);

  % Feedback state    % Desired state
  x = state(1:3);     x_d = state_d(1:3);
  xd = state(4:6);    xd_d = state_d(4:6);
  e = state(7:9);     e_d = state_d(7:9);
  ed = state(10:12);  ed_d = state_d(10:12);

  % Forces and torques, Eqn.[23]
  tau = zeros(3,1);
  Cp = cos(e(1));
  Ct = cos(e(2));

  T =  (g + K(1) * (xd_d(3) - xd(3)) + K(5) * (x_d(3) - x(3))) * m / (Cp * Ct);
  tau(1) = (K(2) * (ed_d(1) - ed(1)) + K(6) * (e_d(1) - e(1))) * Ixx;
  tau(2) = (K(3) * (ed_d(2) - ed(2)) + K(7) * (e_d(2) - e(2))) * Iyy;
  tau(3) = (K(4) * (ed_d(3) - ed(3)) + K(8) * (e_d(3) - e(3))) * Izz;

  % Rotor commands, Eqn.[24]
  w = zeros(3,1);
  w(1) = T / (4 * k) - tau(2) / (2 * k *l) - tau(3) / (4 * b);
  w(2) = T / (4 * k) - tau(1) / (2 * k *l) + tau(3) / (4 * b);
  w(3) = T / (4 * k) + tau(2) / (2 * k *l) - tau(3) / (4 * b);
  w(4) = T / (4 * k) + tau(1) / (2 * k *l) + tau(3) / (4 * b);
  w = sqrt(w);
endfunction

fprintf("Simulation start...\r\n");

% Time params
dt = 0.01;
stop_time = 6;
time = 0:dt:stop_time;

% Physical params
g = 9.81;       % Gravity's acceleration, m/(s*s)
m = 0.468;      % Mass, kg
l = 0.225;      % Rotor to COM distance, m
k = 2.980e-6;   % Lift coefficient
b = 1.140e-7;   % Drag constant
Im = 3.357e-5;  % MOI of rotor, kg*m*m
Ixx = 4.856e-3; % Quadrotor MOI about X-axis, kg*m*m
Iyy = 4.856e-3; % Quadrotor MOI about Y-axis, kg*m*m
Izz = 8.801e-3; % Quadrotor MOI about Z-axis, kg*m*m

% Aerodynamic drag coefficients
Ax = 0.25; % kg/s
Ay = 0.25; % kg/s
Az = 0.25; % kg/s

phy_params = [g, m, l, k, b, Im, Ixx, Iyy, Izz];
A = [Ax, Ay, Az];

% PID gains: kd_z, kd_phi, kd_theta, kd_psi, kp_z, kp_phi, kp_theta, kp_psi
K = [2.5, 1.75, 1.75, 1.75, 1.5, 6, 6, 6];

% Desired state
xi_d = [0, 0, 0];
xi_dot_d = [0, 0, 0];
eta_d = [0, 0, 0];
eta_dot_d = [0, 0, 0];
state_d = [xi_d, xi_dot_d, eta_d, eta_dot_d]';

% Initial conditions
xi = [0, 0, 1]; % m
xi_dot = [0, 0, 0]; % m/s
eta = deg2rad([10, 10, 10]); % rad
eta_dot = [0, 0, 0]; % rad/s

% Memory allocation
state = zeros(12, length(time));
w = zeros(4, length(time));
state(:,1) = [xi, xi_dot, eta, eta_dot]';

% Numerical integration
for t = 1:length(time)-1
  w(:,t) = pd_control(phy_params, state_d, state(:,t), K);
  fn = @(y)quadcoptor(y, phy_params, A, w(:,t));
  state(:,t+1) = rk4(fn, state(:,t), dt);
end
w(:,end) = pd_control(phy_params, state_d, state(:,end), K);

% Control input
figure;
plot(time, w(1,:), 'LineStyle', '-', 'LineWidth', 2); hold on;
plot(time, w(2,:), 'LineStyle', '--', 'LineWidth', 2);
plot(time, w(3,:), 'LineStyle', '-.', 'LineWidth', 2);
plot(time, w(4,:), 'LineStyle', ':', 'LineWidth', 2);
xlabel('Time (s)');
ylabel('Angular velocity (rad/sec)');
title('Control inputs: Rotor angular velocities');
legend('\omega_{1}', '\omega_{2}', '\omega_{3}', '\omega_{4}');
grid on; grid minor;

% Position
figure;
plot(time, state(1,:), 'LineStyle', '-', 'LineWidth', 2); hold on;
plot(time, state(2,:), 'LineStyle', '--', 'LineWidth', 2);
plot(time, state(3,:), 'LineStyle', ':', 'LineWidth', 2);
legend('x', 'y', 'z');
title('Positions x, y, and z');
ylabel('Position (m)');
xlabel('Time (s)');
grid on; grid minor;

% Orientation
figure;
plot(time, rad2deg(state(7,:)), 'LineStyle', '-', 'LineWidth', 2); hold on;
plot(time, rad2deg(state(8,:)), 'LineStyle', '--', 'LineWidth', 2);
plot(time, rad2deg(state(9,:)), 'LineStyle', ':', 'LineWidth', 2);
legend('\phi', '\theta', '\psi');
title('Angles \phi, \theta, and \psi');
xlabel('Time (s)');
ylabel('Angle (degree)');
grid on; grid minor;

fprintf("Simulation complete!\r\n");
