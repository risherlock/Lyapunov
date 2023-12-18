function y_update = rk4(f,y,dt)
  k1 = f(y);
  k2 = f(y + 0.5 * dt * k1);
  k3 = f(y + 0.5 * dt * k2);
  k4 = f(y + dt * k3);

  K = (1 / 6) * (k1 + 2 * k2 + 2 * k3 + k4);
  y_update = y + K * dt;
end