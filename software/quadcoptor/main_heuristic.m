% Numerical simulation of quadrotor dynamics with PD control
% 2023-12-02

clc
clear
close all

clear -global zn_d znd_d

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
endfunction

% Computes roll, pitch, and total thrust from jounce and acceleration.
function w = feedforward_control(phy_params, ddxi, dxi, A, dt)
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

    % Virtual control vector, Eqn.(26)
    dx = ddxi(1) + A(1) * dxi(1);
    dy = ddxi(2) + A(2) * dxi(2);
    dz = ddxi(3) + A(3) * dxi(3);

    psi = 0;
    cy = cos(psi);
    sy = sin(psi);

    phi = asin((dx * sy - cy * dy) / (dx^2 + dy^2 + (dz + g)^2));
    theta = atan2((cy * dx + dy * sy), (dz + g));

    cr = cos(phi);
    sr = sin(phi);
    cp = cos(theta);
    sp = sin(theta);

    % Total thrust, Eqn.(26)
    T = m * (dx * (sp * cy * cr + sy * sr) ...
      + dy * (sp * sy * cr - cy * sr) ...
      + (dz + g) * cp * cr);

    global zn_d = [0, 0, 0]';
    global znd_d = [0, 0, 0]';
    n_d = [phi, theta, 0]';
    nd_d = (n_d - zn_d) / dt;
    ndd_d = (nd_d - znd_d) / dt;
    zn_d = n_d;
    znd_d = nd_d;

    % Eqn.(16)
    W = [1, 0, -sp; 0, cr, cp * sr; 0, -sr, cp * cr];
    J = W' * diag([Ixx, Iyy, Izz]) * W;

    % Eqn.(19)
    Iyz = Iyy - Izz;
    Izy = Izz - Iyy;
    C = zeros(3,3);
    C(1,2) = Iyz * (nd_d(2) * cr * sr + nd_d(3) * sr^2 * cp) ...
           + Izy * nd_d(3) * cr^2 * cp - Ixx * nd_d(3) * cp;
    C(1,3) = Izy * nd_d(3) * cr * sr * cp^2;
    C(2,1) = Izy * (nd_d(2) * cr * sr + nd_d(3) * sr * cp) ...
           + Iyz * nd_d(3) * cr^2 * cp + Ixx * nd_d(3) * cp;
    C(2,2) = Izy * nd_d(1) * cr * sr;
    C(2,3) = -Ixx * nd_d(3) * sp * cp + Iyy * nd_d(3) * sr^2 * sp * cp ...
           + Izz * nd_d(3) * cr^2 * sp * cp;
    C(3,1) = Iyz * nd_d(3) * cp^2 * sr * cr - Ixx * nd_d(2) * cp;
    C(3,2) = Izy * (nd_d(2) * cr * sr * sp + nd_d(1) * sr^2 * cp) ...
           + Iyz * nd_d(1) * cr^2 * cp ...
           + (Ixx - Iyy * sr^2 - Izz * cr^2) * nd_d(3) * sp * cp; 
    C(3,3) = Iyz * nd_d(1) * cr * sr * cp^2 - Iyy * nd_d(2) * sr^2 * cp * sp ...
           - Izz * nd_d(2) * cr^2 * cp * sp + Ixx * nd_d(2) * cp * sp;

   % Torque, Eqn.(20)
   tau = J * ndd_d + C * nd_d;

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
stop_time = 7;
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
  ddxi = [acceleration(t), 0, 0]';
  dxi = [velocity(t), 0, 0]';

  w(:,t) = feedforward_control(phy_params, ddxi, dxi, A, dt);
  fn = @(y)quadcoptor(y, phy_params, A, w(:,t));
  state(:,t+1) = rk4(fn, state(:,t), dt);
end
ddxi = [acceleration(end), 0, 0];
dxi = [velocity(end), 0, 0];
w(:,end) = feedforward_control(phy_params, ddxi, dxi, A, dt);

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
