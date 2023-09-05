% Numerical simulation of quadrotor dynamics with open loop control
% Reference: Modelling and control of quadcoptor (2011) by Teppo Luukkonen
% Rishav (2023-07-26)

clc
clear
close all

fprintf("Simulation start...\r\n");

% Time params
dt = 0.001;
stop_time = 2;
time = 0:dt:stop_time;

% Physical params
g = 9.81;         % Gravity's acceleration, m/(s*s)
m = 0.468;        % Mass, kg
l = 0.225;        % Rotor to COM distance, m
k = 2.980e-6;     % Lift coefficient
b = 1.140e-7;     % Drag constant
Im = 3.357e-5;    % MOI of rotor, kg*m*m
Ixx = 4.856e-3;   % Quadrotor MOI about X-axis, kg*m*m
Iyy = 4.856e-3;   % Quadrotor MOI about Y-axis, kg*m*m
Izz = 8.801e-3;   % Quadrotor MOI about Z-axis, kg*m*m

% Aerodynamic drag coefficients
Ax = 0.25; % kg/s
Ay = 0.25; % kg/s
Az = 0.25; % kg/s

phy_params = [g, m, l, k, b, Im, Ixx, Iyy, Izz];
aero_params = [Ax, Ay, Az];

% Initial conditions
xi = [0, 0, 0];
xi_dot = [0, 0, 0];
eta = [0, 0, 0];
eta_dot = [0, 0, 0];

% Rotor command profile params
t1s = 0;
t2s = 0.5;
t3s = 1;
t4s = 1.5;
t4e = 2;

amp_t1 = 65;
amp_t2 = 35;
amp_t3 = 0;
amp_t4 = -35;
offset = 620;

% Generate rotor command profile
w = zeros(4, length(time));
for i = 1:length(time)
  t = time(i);
  if t >= t1s && t < t2s
    w(1, i) = amp_t1 * sin(2 * pi * (t - t1s) / (t2s - t1s)) + offset;
    w(2, i) = w(1, i);
    w(3, i) = w(1, i);
    w(4, i) = w(1, i);
  elseif t >= t2s && t <= t3s
    w(1, i) = offset;
    w(2, i) = amp_t4 * sin(2 * pi * (t - t2s) / (t3s - t2s)) + offset;
    w(3, i) = offset;
    w(4, i) = amp_t2 * sin(2 * pi * (t - t2s) / (t3s - t2s)) + offset;
  elseif t >= t3s && t <= t4s
    w(1, i) = amp_t4 * sin(2 * pi * (t - t3s) / (t4s - t3s)) + offset;
    w(2, i) = offset;
    w(3, i) = amp_t2 * sin(2 * pi * (t - t3s) / (t4s - t3s)) + offset;
    w(4, i) = offset;
  elseif t >= t4s && t <= t4e
    w(1, i) = amp_t2 * sin(2 * pi * (t - t4s) / (t4e - t4s)) + offset;
    w(2, i) = amp_t4 * sin(2 * pi * (t - t4s) / (t4e - t4s)) + offset;
    w(3, i) = w(1, i);
    w(4, i) = w(2, i);
  end
end

% Memory allocation
state = zeros(12, length(time));
state(:,1) = [xi, xi_dot, eta, eta_dot]';

% Numerical integration
for t = 1:length(time)-1
  fn = @(y)quadcoptor(y, phy_params, aero_params, w(:,t));
  state(:,t+1) = RK4(fn, state(:,t), dt);
end

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
