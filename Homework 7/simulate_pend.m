clear all
close all
clc

%%

c1 = 0.5;
l1 = 1;
I1 = 0.05;
m1 = 1;

c2 = 0.25;
l2 = 0.5;
I2 = 0.05;
m2 = 1;

g = 9.81;

param = [c1; l1; c2; l2; m1; I1; m2; I2; g];

p0 = [0, 0];

z = [pi/6; pi/6; 0; 0];
dim = length(z);
num_step = 4000;
dt = 0.001;
energy_trj = [];

x0 = [0.5; -1];
radi = 0.25;
omega = 2*pi*0.5;

tspan = zeros(1, num_step);
x_des_trj = zeros(dim/2, num_step);
dx_des_trj = zeros(dim/2, num_step);

x_trj = zeros(dim/2, num_step);
dx_trj = zeros(dim/2, num_step);

z_trj = zeros(dim, num_step);

for i=1:num_step
    t = i*dt;
    tspan(i) = t;
   
   x_des = x0 + [radi*cos(omega*t); radi*sin(omega*t)];
   dx_des = [-radi*omega*sin(omega*t); radi*omega*cos(omega*t)];
   ddx_des = [-radi*omega^2*cos(omega*t); -radi*omega^2*sin(omega*t)];
   
   x_des_trj(:,i) = x_des;
   dx_des_trj(:,i) = dx_des;
   
   u = controller(z, param, x_des, dx_des, ddx_des);
   dz = dyn_pend(z, u, param);
   
   
   z(dim/2+1:end) = z(dim/2+1:end) + dz(dim/2+1:end) * dt;
   z(1:dim/2) = z(1:dim/2) + z(dim/2+1:end)*dt;
   z_trj(:,i) = z;
   
   key_pt = keypoints_pend(z, param);
   rA = key_pt(:,1);
   rB = key_pt(:,2);
    
   x_trj(:,i) = rB;
   dx_trj(:,i) = velocity_rB(z, param);
end

figure
subplot(2,1,1)
plot(tspan, x_des_trj(1,:), tspan, x_trj(1,:));
ylabel('x')
subplot(2,1,2)
plot(tspan, x_des_trj(2,:), tspan, x_trj(2,:));
ylabel('y')

%% 
figure
% Prepare plot handles
hold on

TH = 0:.1:2*pi;
plot( x0(1) + radi * cos(TH), ...
      x0(2) + radi * sin(TH),'k--'); 
% plot(x0(1), x0(2),'*')
h_OA = plot([0],[0],'LineWidth',4);
h_AB = plot([0],[0],'LineWidth',4);


xlabel('x'); ylabel('y');
h_title = title('t=0.0s');

axis equal
axis([-1.5 1.5 -1.5 1.5]);

%Step through and update animation
for i = 1:length(tspan)
    % skip frame.
    if mod(i,10)
        continue;
    end
    t = tspan(i);
    z = z_trj(:,i); 
    keypoints = keypoints_pend(z,param);

    rA = keypoints(:,1); % Vector to base of cart
    rB = keypoints(:,2);

    set(h_title,'String',  sprintf('t=%.2f',t) ); % update title

    set(h_OA,'XData',[0 rA(1)]);
    set(h_OA,'YData',[0 rA(2)]);

    set(h_AB,'XData',[rA(1) rB(1)]);
    set(h_AB,'YData',[rA(2) rB(2)]);

    pause(.01)
end
