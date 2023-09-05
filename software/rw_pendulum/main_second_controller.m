% Rishav (2021-08-06)

clc
clear
close all

% ~~ Configuration start ~~%

% Physical parameters
m1 = 0.02;     % Mass of pendulum, kg
m2 = 0.3;      % Mass of wheel, kg
l1 = 0.123;    % Length of the pendulum, m
lc1 = 0.063;   % Pivot to COM distance, m
I1 = 47e-6;    % MOI of pendulum, Kg*m^2
I2 = 32e-6;    % MOI of wheel, Kg*m^2
tau = 0;       % Torque input applied to the disk
g = 9.804;     % Acceleration due to gravity, m/s^2

% Initial Conditions
init_q1 = 0.8*pi;    % Angle between vertical and pendulum, rad
init_q2 = 0;         % Angle of the wheel, rad

init_q1_dot = 0;     % rad/s
init_q2_dot = 0;     % rad/s

% Simulation params
dt = 0.001;         % Sample time, s
start_time = 0;     % s
stop_time = 300;    % s

% Control gains
ke = 400;
kv = 0.01;
kp = 0.1;
kd = 0.05;

% ~~ Configuration end ~~ %

% Inertia matrix: Constant and positive definite
D = [m1 * lc1^2 + m2 * l1^2 + I1 + I2, I2; I2, I2];

% Intermediate variables
m_ = m1 * lc1 + m2 * l1;
k1 = kv * D(2,1) * m_ * g / det(D);
k2 = kv * D(1,1) / det(D);

time = start_time:dt:stop_time;
state = zeros(4, length(time));
state(:,1) = [init_q1, init_q1_dot, init_q2, init_q2_dot];
E_plot = zeros(1, length(time)-1);   % Energy function
V_plot = zeros(1, length(time)-1);   % Lyapunov function
tau_plot = zeros(1, length(time)-1); % Command torque to motor

disp("Simulation in progress...");
for i_iters = 1: length(time) - 1
  q_dot = [state(2,i_iters), state(4,i_iters)]';
  q = [state(1,i_iters), state(3,i_iters)]';
    
  % Total energy of pendulum
  E = 0.5 * q_dot' * D * q_dot + m_ * g * (cos(state(1,i_iters)) - 1); 
    
  % Reaction wheel control torque
  tau = (- kd * q_dot(2) - kp * q(2) + k1 * sin(q(1))) / (ke * E + k2); 
    
  % Integrate the dynamics of RWP
  fn = @(y)rwp_dynamics(y, m_ * g, D, tau);
  state(:,i_iters+1) = rk4(fn, state(:,i_iters), dt);
    
  V_plot(i_iters) = 0.5 * ke * E^2 + ...
                    0.5 * kv * q_dot(2)^2 + ...
                    0.5 * kp * q(2)^2;
  E_plot(i_iters) = E;
  tau_plot(i_iters) = tau;
end

disp("Simulation complete !")
plot_output(E_plot, tau_plot, V_plot, state, time);
