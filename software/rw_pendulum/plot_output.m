function [] = plot_output(E_plot, tau_plot, V_plot, state, time)
  % Energy, Command, Lapynov fn and phase plot
  subplot(2,2,1);
  plot(time(1:length(time)-1), E_plot, "lineWidth", 1);
  title("Energy function");
  xlabel("Time (s)");
  grid on;
    
  subplot(2,2,2);
  plot(time(1:length(time)-1), tau_plot, "lineWidth", 1);
  xlabel("Time (s)");
  title("Control input");
  grid on;
    
  subplot(2,2,3);
  plot(time(1:length(time)-1), V_plot, "lineWidth", 1);
  title("Lyapunov function");
  xlabel("Time (s)");
  grid on;
    
  subplot(2,2,4);
  plot(state(1,:), state(2,:), "lineWidth", 1);
  title("Plase plot");
  xlabel("dq_1")
  ylabel("q_1");
  grid on;
    
  % q1, dq1, q2 and dq2
  figure;
  subplot(2,2,1);
  plot(time, state(1,:), "lineWidth", 1);
  title("Angle of the pendulum");
  ylabel("q_1 (rad)")
  xlabel("Time (s)");
  grid on;
    
  subplot(2,2,2);
  plot(time, state(2,:),"lineWidth", 1);
  title("Angular velocity of the pendulum");
  ylabel("q_1dot (rad/s)")
  xlabel("Time (s)");
  grid on;
    
  subplot(2,2,3);
  plot(time, state(3,:), "lineWidth", 1);
  title("Angle of the reaction wheel");
  ylabel("q_2 (rad)")
  xlabel("Time (s)");
  grid on;
    
  subplot(2,2,4);
  plot(time, state(4,:), "lineWidth", 1);
  title("Anglular velocity of the reaction wheel");
  ylabel("q_2dot (rad/s)")
  xlabel("Time (s)");
  grid on;
end
