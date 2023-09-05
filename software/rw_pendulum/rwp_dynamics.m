function [x_dot] = rwp_dynamics(x, m_xg, D, tau)
%%% Differential equation governing reaction wheel pendulum
%   
% Inputs:
%   x       = State vector: [q1, q1_dot, q2, q2_dot]'
%   m_xg    = (m1 * lc1 + m2 * l1) * g
%   D       = [m1 * lc1^2 + m2 * l1^2 + I1 + I2, I2; I2, I2]
%   tau     = Reaction wheel control torque
%
% Output:
%   x_dot = Derivative of state vector
%
% Reference: Isabelle Fantoni et al. 
%   [1] Stabilization of the reaction wheel pendulum using an energy approach (2001)
%   [2] Non-linear control for underactuated mechanical systems - Chapter 7 (2002)

  q1 = x(1);
  q1_dot = x(2);
  q2_dot = x(4);
    
  % RWP dynamics
  u = [0, tau]'; % Joint torques
  g_q = [-m_xg * sin(q1), 0]';
  q_dot_dot = D\(u - g_q);

  x_dot = [q1_dot, q_dot_dot(1), q2_dot, q_dot_dot(2)]';
end
