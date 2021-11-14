function u = controller(z, param, x_des, dx_des, ddx_des) 
   % ******** Implement your controller ********
   ddt = @(r) jacobian(r, [q; dq])*[dq; ddq];
   dim = length(z);
   key_pt = keypoints_pend(z, param);
   rB = key_pt(1:dim/2, 2);
   rB_dot = velocity_rB(z, param);
   
   position_err = x_des - rB;
   velocity_err = dx_des - rB_dot;
   
   K = [50, 0; 0, 50];
   D = [5, 0; 0, 5];
   
   command = ddx_des + K * position_err + D * velocity_err;
   command2 = K * position_err + D * velocity_err;
   
   u = zeros(2,1); % torque
   J = Jacobian_rB(z, param); %jacobian function
   J_dot = Jdot_rB(z, param);
   M = A_pend(z, param); % mass matrix
   z_zero_vel = z;
   z_zero_vel(dim/2 + 1:end) = zeros(dim/2, 1);
   u_zero = zeros(size(command));
   G = -b_pend(z_zero_vel, u_zero, param); % gravitational part
   C = -b_pend(z, u_zero, param) - G; % coriolis part
   
   lambda_inv = J * inv(M) * J.';
   lambda = inv(lambda_inv);
   mu = lambda * J * inv(M) * C - lambda * J_dot * z(dim/2 + 1:end);
   rho = lambda * J * inv(M) * G;
   
%    F = lambda * command + mu + rho;
%    F = lambda * command + mu;
   F = lambda * command + rho;
%    F = lambda * command2 + mu + rho;
   
   u = J.' * F;
   
end