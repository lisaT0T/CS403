clear all
close all
clc

syms th1 th2 dth1 dth2 ddth1 ddth2 real
syms c1 c2 m1 m2 l1 l2 i1 i2 tau1 tau2 g real

q = [th1; th2];
dq = [dth1; dth2];
ddq = [ddth1; ddth2];

u = [tau1; tau2];
p = [c1 c2 m1 m2 l1 l2 i1 i2 g];

rc1 = [c1*sin(th1); -c1*cos(th1)]; % position at center of mass link1
rB = [l1*sin(th1); -l1*cos(th1)]; % position at point B
rc2 = [l1*sin(th1) + c2*sin(th1+th2); -l1*cos(th1) - c2*cos(th1+th2)]; % position at center of mass link2
rC = [l1*sin(th1) + l2*sin(th1+th2); -l1*cos(th1) - l2*cos(th1+th2)]; % end effector (position at point A)

ddt = @(r) jacobian(r, [q; dq])*[dq; ddq]; % function to derive the velocity and acceleration

vc1 = ddt(rc1);
vc2 = ddt(rc2);

%kinetic energy
T = 1/2 * m1 * dot(vc1, vc1) + 1/2 * m2 * dot(vc2, vc2) + 1/2 * i1 * dth1^2 + 1/2 * i2 * dth2^2; 

h1 = dot(rc1, -[0; -1]);
h2 = dot(rc2, -[0; -1]);

V = m1*g*h1 + m2*g*h2; %potential energy

L = T - V; % lagrangian
E = T + V; % total energy

Q_tau = [tau1; tau2];
g1 = ddt(jacobian(L, dq)') - jacobian(L, q)' - Q_tau;

A = simplify(jacobian(g1, ddq));
b = simplify(A*ddq - g1);

z = [q; dq];

matlabFunction(A,'file','A_pend','vars',{z p});
matlabFunction(b,'file','b_pend','vars',{z u p});
matlabFunction(rB,'file','point_b','vars',{z p});
matlabFunction(rC,'file','end_pend','vars',{z p});
matlabFunction(E,'file','energy_pend','vars',{z p});

% simulation
th1 = 3;
th2 = 0;
dth1 =0;
dth2 = 0;
m1 = 1;
m2 = 1;
i1 = 0.05;
i2 = 0.05;
l1 = 1;
l2 = 0.5;
c1 = 0.5;
c2 = 0.25;
g = 9.81;
tau1 = 0;
tau2 = 0;

q = [th1; th2];
dq = [dth1; dth2];

% u = [tau1; tau2];
z = [q; dq];
param = [c1 c2 m1 m2 l1 l2 i1 i2 g];

p0 = [0; 0];

rB = [l1*sin(th1); -l1*cos(th1)];
rC = [l1*sin(th1) + l2*sin(th1+th2); -l1*cos(th1) - l2*cos(th1+th2)];

dt = 0.001;
E_trj = [];
th1_graph = [];
th2_graph = [];
for i = 1:7000 % time interval [0s, 7s]
   u = [0; 0];
   dz = dynamics(z, u, param); % update on [dq, ddq]
   
   z(3:4) = dz(1:2) + dz(3:4)*dt; % update dq
   z(1:2) = z(1:2) + z(3:4)*dt; % update

   E_trj(i) = energy_pend(z, param);
   th1_graph(i) = z(1);
   th2_graph(i) = z(2);
   
   pB = point_b(z, param); 
   pC = end_pend(z, param);
   
   if mod(i, 10) == 1
       drawLine2D(p0, pB);
       hold on
       drawLine2D(pB, pC);
       hold  off
       axis equal
       xlim([-1.5, 1.5]);
       ylim([-1.5, 1.5]);
       pause(0.0001);
   end

end

figure
plot(E_trj);
ylim([15, 20]);

figure
plot(th1_graph);

figure
plot(th2_graph);

function dz = dynamics(z, u, param)
    A = A_pend(z, param);
    b = b_pend(z, u, param);
    
    ddq = inv(A)*b;
    
    % dz = [dq, ddq]
    dz = zeros(4,1);
    dz(1:2) = z(3:4);
    dz(3:4) = ddq;    
end

function x = drawLine2D(start_pt, end_pt)
    % draw the line starting from 'start_pt' and ending at 'end_pt'
    % start_pt = [x;y]
    % end_pt = [x;y]

    plot([start_pt(1), end_pt(1)], [start_pt(2), end_pt(2)], 'linewidth',4);
end