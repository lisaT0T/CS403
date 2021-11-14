clear all

%% Paramter preparation
syms th1 dth1 ddth1 th2 dth2 ddth2 real
syms c1 l1 c2 l2 m1 I1 m2 I2 real
syms g tau1 tau2 real

q = [th1; th2];
dq = [dth1; dth2];
ddq = [ddth1; ddth2];

u = [tau1; tau2];
p = [c1; l1; c2; l2; m1; I1; m2; I2; g];

%% Kinetic and Potential engergy
Rc1 = [c1*sin(th1); -c1*cos(th1)];
Rc2 = [l1*sin(th1) + c2*sin(th1+th2); -l1*cos(th1)-c2*cos(th1+th2)];

ddt = @(r) jacobian(r, [q; dq])*[dq; ddq];

% ******** Define Kinetic and Potential Energy (T, V) ********
Vc1 = ddt(Rc1);
Vc2 = ddt(Rc2); %r dot(velocity) of c1 c2

T = 0.5 * m1 * dot(Vc1, Vc2) + 0.5 * m2 * dot(Vc2, Vc2) + 0.5 * I1 * dth1^2 + 0.5 * I2 * dth2^2;

h1 = dot(Rc1, -[0; -1]);
h2 = dot(Rc2, -[0; -1]); % h = -g
V = m1 * g * h1 + m2 * g * h2;

%% Generalized Force

% ******** Define Generalized Force (Q) ********
w = [dth1; dth2];
J = jacobian(w, dq);

Q = J.'*u; % Q = [tau1; tau2]

%% Lagrange equation
L = T - V;
g = ddt(jacobian(L,dq).') - jacobian(L,q).' - Q;

A = simplify(jacobian(g,ddq));
b = simplify(A*ddq - g);

gravity = simplify(jacobian(V,q)).';
coriolis = simplify(-b - gravity + Q);

z = [q; dq];
dz = [dq; ddq];

rA = [l1*sin(th1); -l1*cos(th1)];
rB = [l1*sin(th1) + l2*sin(th1+th2); -l1*cos(th1) - l2*cos(th1+th2)];
keypoints = [rA rB];

J_B = jacobian(rB, q);
J_B_dot = reshape(ddt(J_B(:)), size(J_B));
vB = J_B*dq;

%% Save files
matlabFunction(A, 'file', 'A_pend', 'var', {z p});
matlabFunction(b, 'file', 'b_pend', 'var', {z u p});
matlabFunction(keypoints, 'file','keypoints_pend', 'var',{z p});
matlabFunction(J_B, 'file', 'Jacobian_rB','var',{z p});
matlabFunction(vB, 'file', 'velocity_rB','var',{z p});
matlabFunction(gravity, 'file','grav_pend', 'var', {z p});
matlabFunction(coriolis, 'file','coriolis_pend', 'var', {z p});
matlabFunction(J_B_dot, 'file', 'Jdot_rB', 'var',{z p});