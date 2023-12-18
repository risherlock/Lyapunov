% Numerical simulation of quadrotor dynamics with PD control
% 2023-12-08

clc
clear
close all

clear -global zn_d zxd

% Simple integration of time series
function x = integrate(dx, dt)
  x = zeros(size(dx));
  x(1) = dx(1);
  for i = 2:length(x)
    x(i) += x(i-1) + dx(i);
  end
  x = x .* dt;
endfunction

% Heuristic jounce generator Eqn.(28)
function j = get_heuristic_jounce(a, b, c, time)
  j = zeros(size(time));
  sign_flag = 1;

  for i = 1 : length(time)
    t = time(i);
    if(t <= 7)
      while(t > 2 * c)
        t = t - 2 * c;
        sign_flag = 0;
      end

      factor = pi * t / b;
      if (t <= b)
        j(i) = a * sin(factor);
      elseif (t <= 3 * b)
       j(i) = -a * sin(0.5 * (factor - pi));
      elseif (t <= 4 * b)
        j(i) = a * sin(factor - 3 * pi);
      end

      if(!sign_flag)
        j(i) = -j(i);
      end
    end
  end
end


% Computes roll, pitch, and total thrust from jounce and acceleration.
function w = feedforward_control(phy_params, state, xdd_d, xd_d, x_d, K, dt)
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

    % Feedback state
    x = state(1:3);
    xd = state(4:6);
    n = state(7:9);
    nd = state(10:12);

    % Virtual control vector, Eqn.(26)
    global zxd = [0, 0, 0]';
    xdd = (xd - zxd) / dt;
    zxd = xd;

    dx = K(1) * (x_d(1) - x(1)) + K(7) * (xd_d(1) - xd(1)) + K(13) * (xdd_d(1) - xdd(1));
    dy = K(2) * (x_d(1) - x(2)) + K(8) * (xd_d(1) - xd(2)) + K(14) * (xdd_d(1) - xdd(2));
    dz = K(3) * (x_d(1) - x(3)) + K(9) * (xd_d(1) - xd(3)) + K(15) * (xdd_d(1) - xdd(3));
    psi_d = 0;
    cy = cos(psi_d);
    sy = sin(psi_d);

    phi_d = asin((dx * sy - cy * dy) / (dx^2 + dy^2 + (dz + g)^2));
    theta_d = atan((cy * dx + dy * sy) / (dz + g));

    cr = cos(phi_d);
    sr = sin(phi_d);
    cp = cos(theta_d);
    sp = sin(theta_d);

    % Total thrust, Eqn.(26)
    T = m * (dx * (sp * cy * cr + sy * sr) ...
      + dy * (sp * sy * cr - cy * sr) ...
      + (dz + g) * cp * cr);

    global zn_d = [0, 0, 0]';
    n_d = [phi_d, theta_d, psi_d]';
    nd_d = (n_d - zn_d) / dt;
    zn_d = n_d;

    % Torques, Eqn.[30]
    tau = zeros(3,1);
    tau(1) = (K(4) * (n_d(1) - n(1)) + K(10) * (nd_d(1) - nd(1))) * Ixx;
    tau(2) = (K(5) * (n_d(2) - n(2)) + K(11) * (nd_d(2) - nd(2))) * Iyy;
    tau(3) = (K(6) * (n_d(3) - n(3)) + K(12) * (nd_d(3) - nd(3))) * Izz;

    % Rotor angular velocity
    w = zeros(4,1);
    w(1) = T / (4 * k) - tau(2) / (2 * k * l) - tau(3) / (4 * b);
    w(2) = T / (4 * k) - tau(1) / (2 * k * l) + tau(3) / (4 * b);
    w(3) = T / (4 * k) + tau(2) / (2 * k * l) - tau(3) / (4 * b);
    w(4) = T / (4 * k) + tau(1) / (2 * k * l) + tau(3) / (4 * b);
    w = sqrt(w);
endfunction
  
fprintf("Simulation start...\r\n");

% Time params
dt = 0.01;
stop_time = 14;
time = 0:dt:stop_time;

% Physical params
g = 9.81;       % Gravity's acceleration, m/s^2
m = 0.468;      % Mass, kg
l = 0.225;      % Rotor to COM distance, m
k = 2.980e-6;   % Lift coefficient, N*s^2 / rad^2
b = 1.140e-7;   % Drag constant, Nm*s^2 / rad^2
Im = 3.357e-5;  % MOI of rotor, kg*m*m
Ixx = 4.856e-3; % Quadrotor MOI about X-axis, kg*m^2
Iyy = 4.856e-3; % Quadrotor MOI about Y-axis, kg*m^2
Izz = 8.801e-3; % Quadrotor MOI about Z-axis, kg*m^2

% Aerodynamic drag coefficients
Ax = 0.25; % kg/s
Ay = 0.25; % kg/s
Az = 0.25; % kg/s

phy_params = [g, m, l, k, b, Im, Ixx, Iyy, Izz];
A = [Ax, Ay, Az];

%  1-6 : kp_x, kp_y, kp_z, kp_phi, kp_theta, kp_psi
%  7-12: kd_x, kd_y, kd_z, kd_phi, kd_theta, kd_psi
% 13-15: kdd_x, kdd_y, kdd_z
K = [1.85, 8.55, 1.85, 3, 3, 3, 0.75, 0.75, 0.75, 0.75, 0.75, 0.75, 1, 1, 1];

% Trajectory generation
a = 1;
b = 0.5;
c = 2;
jounce = get_heuristic_jounce(a, b, c, time);
jerk = integrate(jounce, dt);
acceleration = integrate(jerk, dt);
velocity = integrate(acceleration, dt);
position = integrate(velocity, dt);

% Initial conditions
xi = [0, 0, 0]; % m
xi_dot = [0, 0, 0]; % m/s
eta = deg2rad([0, 0, 0]); % rad
eta_dot = [0, 0, 0]; % rad/s

% Memory allocation
state = zeros(12, length(time));
w = zeros(4, length(time));
state(:,1) = [xi, xi_dot, eta, eta_dot]';

% Numerical integration
for t = 1:length(time)-1
  xdd_d = [acceleration(t), 0, 0]';
  xd_d = [velocity(t), 0, 0]';
  x_d = [position(t), 0, 0]';
  
  w(:,t) = feedforward_control(phy_params, state(:,t), xdd_d, xd_d, x_d, K, dt);
  fn = @(y)quadcoptor(y, phy_params, A, w(:,t));
  state(:,t+1) = rk4(fn, state(:,t), dt);
  
  temp = quadcoptor(state(:,t), phy_params, A, w(:,t));
  xdd = temp(4:6); 
end
xdd_d = [acceleration(end), 0, 0]';
xd_d = [velocity(end), 0, 0]';
x_d = [position(end), 0, 0]';
w(:,end) = feedforward_control(phy_params, state(:,t), xdd_d, xd_d, x_d, K, dt);

figure;
plot(time, jounce, 'LineWidth', 2);
title("Jounce");
grid on; grid minor;

figure;
plot(time, jerk, 'LineWidth', 2);
title("Jerk");
grid on; grid minor;

figure;
plot(time, acceleration, 'LineWidth', 2);
title("Acceleration");
grid on; grid minor;

figure;
plot(time, velocity, 'LineWidth', 2);
title("Velocity");
grid on; grid minor;

figure;
plot(time, position, 'LineWidth', 2);
title("Position");
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

fprintf("Simulation complete!\r\n");
