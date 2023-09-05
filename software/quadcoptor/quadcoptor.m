% Dynamical equations of quadcoptor
function [state_dot] = quadcoptor(state, phy_params, aero_params, w)
  % Unpack state variables
  xi = state(1:3);
  xi_dot = state(4:6);
  n = state(7:9);
  n_dot = state(10:12);
  
  % Unpack physical params
  g = phy_params(1);
  m = phy_params(2);
  l = phy_params(3);
  k = phy_params(4);
  b = phy_params(5);
  Im = phy_params(6);
  Ixx = phy_params(7);
  Iyy = phy_params(8);
  Izz = phy_params(9);
  
  % Unpack aero-drag params
  Ax = aero_params(1);
  Ay = aero_params(2);
  Az = aero_params(3);
  
  % Angle transedental functions
  cr = cos(n(1)); % roll
  sr = sin(n(1));
  cp = cos(n(2)); % pitch
  sp = sin(n(2));
  tp = tan(n(2));
  cy = cos(n(3)); % yaw
  sy = sin(n(3));
  
  % *** Translational dynamics *** %
  
  % Body Z-axis thrust from rotor angular speed
  T = k*(w(1)^2 + w(2)^2 + w(3)^2 + w(4)^2);
  
  % Forces in body frame acting on quadcoptor
  f_gravity = -g * m * [0, 0, 1]';
  f_thrust = T * [cy*sp*cr + sy*sr; sy*sp*cr - cy*sr; cp*cr];
  f_aero = -diag([Ax, Ay, Az]) * xi_dot;
  
  % Double derivative of inertial position
  xi_ddot = (f_gravity + f_thrust + f_aero) / m;
  
  % *** Rotational dynamics *** %
  
  % Transformation matrix, its inverse, and derivative of the inverse
  W = [1, 0, -sp; 0, cr, cp * sr; 0, -sr, cp * cr];
  Winv = [1, sr * tp, cr * tp; 0, cr, -sr; 0, sr / cp, cr / cp];
  dWinv = zeros(3, 3);
  dWinv(1,2) = n_dot(1)*cr*tp + n_dot(2)*sr/cp^2;
  dWinv(1,3) = -n_dot(1)*sr*cp + n_dot(2)*cr/cp^2;
  dWinv(2,2) = -n_dot(1)*sr;
  dWinv(2,3) = -n_dot(1)*cr;
  dWinv(3,2) = n_dot(1)*(cr  + sr*tp)/cp;
  dWinv(3,3) = (-n_dot(1)*sr + n_dot(2)*cr*tp)/cp;
  
  % Angular velocity in body frame
  v = W * n_dot;
  
  % External and gyro torques in body frame from w
  tau = zeros(3, 1);
  tau(1) = l * k * (w(4)^2 - w(2)^2 );
  tau(2) = l * k * (w(3)^2 - w(1)^2);
  tau(3) = b * (-w(1)^2 + w(2)^2 - w(3)^2 + w(4)^2);
  tau_gyro = Im * (w(1) - w(2) + w(3) - w(4)) * [v(2), -v(1), 0]';
  
  % Angular acceleration in body frame
  I = diag([Ixx, Iyy, Izz]);
  v_dot = I \ (-cross(v, I * v) + tau - tau_gyro);
    
  % Double derivative of inertial angle
  n_ddot = dWinv * v + Winv * v_dot;
  
  % Differential of the state vector
  state_dot = [xi_dot; xi_ddot; n_dot; n_ddot];
end
