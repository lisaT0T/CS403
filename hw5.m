clear all
close all
clc

syms th1 th2 th3 th4 th5 th6;

% each joint consider as a frame
% the position of each frame respect to its last frame
p = {};
p{1} = [0; 0; 0.15];
p{2} = [0.3; 0; 0];
p{3} = [0.15; 0; 0];
p{4} = [0.1; 0; 0];
p{5} = [0.07; 0; 0];
p{6} = [0.05; 0; 0];

S = {};
S{1} = zeros(6, 1); S{1}(3) = 1;
S{2} = zeros(6, 1); S{2}(2) = 1;
S{3} = zeros(6, 1); S{3}(2) = 1;
S{4} = zeros(6, 1); S{4}(1) = 1;
S{5} = zeros(6, 1); S{5}(2) = 1;
S{6} = zeros(6, 1); S{6}(1) = 1;

% desired joint position
T0des = [0 -1 0 0.2;
        1 0 0 0.31;
        0 0 1 0.2;
        0 0 0 1];

th1 = 0;
th2 = 0;
th3 = 0;
th4 = 0;
th5 = 0;
th6 = 0;

% th1 = 0.8;
% th2 = 0.4;
% th3 = -pi/2 - 0.4;
% th4 = -0.9;
% th5 = pi/2 + 0.3;
% th6 = 0.7;
    
th = [th1; th2; th3; th4; th5; th6];

figure
num_step = 50;

for step = 1:num_step
    clf
    T0 = Tset(th, p);
    err = zeros(6, 1); % error between current T0ee with T0desired
    errSO3 = T0des(1:3, 1:3) * transpose(T0{7}(1:3, 1:3));
    so3_skew = logm(errSO3);
    err(1:3) = [so3_skew(3, 2); so3_skew(1, 3); so3_skew(2, 1)];
    err(4:6) = T0des(1:3, 4) - T0{7}(1:3, 4);
    
    J = jacobian(th, p, S);
    th = th + 0.2 * pinv(J) * err; % updata on theta
    
    if norm(err) < 0.001
        break
    end
    
    simulation(th, p);
    drawCoordinate3DScale(T0des(1:3, 1:3), T0des(1:3, 4), 0.1);
    hold on
    grid on
    axis equal  
    view(50, 30);
    hold off
    
    pause(0.001);
end



function J = jacobian(th, p, S)
    % th is an array of theta for each joints
    % p is a set of 
    T0 = Tset(th, p);
    % T0 represent a set of SE(3) wrt to frame 0 
    eeJee = zeros(6, 6); % jacobian in local frame
    for i = 1:6
        Teei = seInv(T0{7}) * T0{i};
        R = Teei(1:3, 1:3);
        p = Teei(1:3, 4);
        p_skew = [0, -p(3), p(2);
                  p(3), 0, -p(1);
                  -p(2), p(1) 0];
        adjoint = zeros(6, 6);
        adjoint(1:3, 1:3) = R;
        adjoint(4:6, 4:6) = R;
        adjoint(4:6, 1:3) = p_skew * R;
        eeJee(:, i) = adjoint * S{i};
    end
    R0ee = T0{7}(1:3, 1:3);
    adjs = zeros(6, 6);
    adjs(1:3, 1:3) = R0ee;
    adjs(4:6, 4:6) = R0ee;
    J = adjs * eeJee; % jacobian in global frame
end

function x = simulation(th, p)
    T0 = Tset(th, p);
    drawLine3D(zeros(1, 3), T0{1}(1:3, 4));
    drawCoordinate3DScale(T0{1}(1:3, 1:3), T0{1}(1:3, 4), 0.05);
    drawLine3D(T0{1}(1:3, 4), T0{2}(1:3, 4));
    drawCoordinate3DScale(T0{2}(1:3, 1:3), T0{2}(1:3, 4), 0.05);
    drawLine3D(T0{2}(1:3, 4), T0{3}(1:3, 4));
    drawCoordinate3DScale(T0{3}(1:3, 1:3), T0{3}(1:3, 4), 0.05);
    drawLine3D(T0{3}(1:3, 4), T0{4}(1:3, 4));
    drawCoordinate3DScale(T0{4}(1:3, 1:3), T0{4}(1:3, 4), 0.05);
    drawLine3D(T0{4}(1:3, 4), T0{5}(1:3, 4));
    drawCoordinate3DScale(T0{5}(1:3, 1:3), T0{5}(1:3, 4), 0.05);
    drawLine3D(T0{5}(1:3, 4), T0{6}(1:3, 4));
    drawCoordinate3DScale(T0{6}(1:3, 1:3), T0{6}(1:3, 4), 0.05);
    drawLine3D(T0{6}(1:3, 4), T0{7}(1:3, 4));
    drawCoordinate3DScale(T0{7}(1:3, 1:3), T0{7}(1:3, 4), 0.05);
    
    text(0, 0, 0, 'frame{0}');
end

function T0 = Tset(th, p)
    R1 = rz(th(1));
    R2 = ry(th(2));
    R3 = ry(th(3));
    R4 = rx(th(4));
    R5 = ry(th(5));
    R6 = rx(th(6));

    % SE3 for each frame wrt its last frame
    T01 = SE3(R1, [0; 0; 0]);
    T12 = SE3(R2, p{1});
    T23 = SE3(R3, p{2});
    T34 = SE3(R4, p{3});
    T45 = SE3(R5, p{4});
    T56 = SE3(R6, p{5});
    T6ee = SE3(eye(3), p{6});

    % SE3 wrt frame 0
    T0 = {};
    T0{1} = T01;
    T0{2} = T01*T12;
    T0{3} = T0{2} * T23;
    T0{4} = T0{3} * T34;
    T0{5} = T0{4} * T45;
    T0{6} = T0{5} * T56;
    T0{7} = T0{6} * T6ee; % end effector
end

function SE = SE3(R, p)
    SE(1:3, 1:3) = R;
    SE(1:3, 4) =  p;
    SE(4, :) = [0 0 0 1];
end

function Rx = rx(theta)
    Rx = [1 0 0; 
          0 cos(theta) -sin(theta); 
          0 sin(theta) cos(theta)];
end

function Ry = ry(theta)
    Ry = [cos(theta) 0 sin(theta); 
          0 1 0; 
          -sin(theta) 0 cos(theta)];
end

function Rz = rz(theta)
    Rz = [cos(theta) -sin(theta) 0; 
          sin(theta) cos(theta) 0;
          0 0 1];
end

function x = drawCoordinate3DScale(R, p, scale)
    hold on
    R = scale*R;
    quiver3(p(1),p(2),p(3), R(1,1), R(2,1), R(3,1),0, 'linewidth',3, 'color','r', 'MaxHeadSize',0.8)
    quiver3(p(1),p(2),p(3), R(1,2), R(2,2), R(3,2),0, 'linewidth',3, 'color','b', 'MaxHeadSize',0.8)
    quiver3(p(1),p(2),p(3), R(1,3), R(2,3), R(3,3),0, 'linewidth',3, 'color','g', 'MaxHeadSize',0.8)
end

function x = drawLine3D(start_pt, end_pt)
    plot3([start_pt(1), end_pt(1)], [start_pt(2), end_pt(2)], [start_pt(3), end_pt(3)], 'linewidth',10);
end

%the inverse function for SE(3)
function inv = seInv(SE)
   %extract the rotation matrix from the given SE(3)
   R = SE(1:3, 1:3);
   %extract p vector from the SE(3)
   p = SE(1:3, 4);
   %transpose the rotation matrix
   Rt = transpose(R);
   RtP = -Rt * p;
   inv(1:3, 1:3) = Rt;
   inv(1:3, 4) = RtP;
   inv(4, :) = [0 0 0 1];
end