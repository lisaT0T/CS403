clear all
close all
clc

addpath('E:\College\2021 Fall\CS 403\classCode\matlab_utils');

% euler angles
eul1 = [0.3 0.2 0.5];
eul2 = [0.7 pi pi/2];
eul3 = [pi/3 0 0];

% Orientation Matrix from the function oriM
R01 = oriM(eul1);
R12 = oriM(eul2);
R23 = oriM(eul3);

% Translation vector
p0 = [0 0 0];
p01 = [0.4 0.8 1.2];
p12 = [-0.4 0.5 1.0];
p23 = [0.5 -0.8 1.2];

% SE3 Matrix from se3 function
T01 = se3(R01, p01);
T12 = se3(R12, p12);
T23 = se3(R23, p23);
T02 = T01*T12;
T03 = T01*T12*T23;

figure
grid on

drawCoordinate3DScale(eye(3), p0, 0.5);
text(0.15, 0.15, 0, 'frame{0}');
drawCoordinate3DScale(T01(1:3, 1:3), T01(1:3, 4), 0.5);
text(T01(1, 4)+0.25, T01(2, 4)+0.25, T01(3, 4)-0.05, 'frame{1}');
drawCoordinate3DScale(T02(1:3, 1:3), T02(1:3, 4), 0.5);
text(T02(1, 4)+0.05, T02(2, 4)+0.05, T02(3, 4), 'frame{2}');
drawCoordinate3DScale(T03(1:3, 1:3), T03(1:3, 4), 0.5);
text(T03(1, 4), T03(2, 4), T03(3, 4)+0.1, 'frame{3}');

xlim([-0.5, 1.5]);
ylim([-1, 1.5]);
zlim([-0.2, 2.5]);

xlabel('x', 'fontsize',22);
ylabel('y', 'fontsize',22);
zlabel('z', 'fontsize',22);
view(60, 40);

figure
t = 0
for i = 1:50
    clf
    t = t + 0.1;
    theta = pi * t;
    %the global frame
    drawCoordinate3DScale(eye(3), zeros(1, 3), 0.5);
    text(0.05, 0.05, 0.4, 'frame{0}');
    text(T03(1, 4), T03(2, 4), T03(3, 4)+0.1, 'frame{3}');
    position = T03(1:3, 4) + [0.1*sin(theta)+0.05, 0.3*cos(theta)+0.08, sin(theta)+0.5];
    drawCoordinate3DScale(T03(1:3, 1:3), position, 0.5);
    
    hold on
    grid on
    
    xlim([-0.5, 1.75]);
    ylim([-1, 1]);
    zlim([-0.2, 2.75]);
    view(170, 15);
    
    xlabel('x', 'fontsize',22);
    ylabel('y', 'fontsize',22);
    zlabel('z', 'fontsize',22);
    pause(.01);
    hold off
end

function R = oriM(x)
    Rz = [cos(x(1)) -sin(x(1)) 0; 
          sin(x(1)) cos(x(1)) 0;
          0 0 1];
    Ry = [cos(x(2)) 0 sin(x(2)); 
          0 1 0; 
          -sin(x(2)) 0 cos(x(2))];
    Rx = [1 0 0; 
          0 cos(x(3)) -sin(x(3)); 
          0 sin(x(3)) cos(x(3))];
    R = Rz * Ry * Rx;
end

function SE = se3(R, p)
    SE(1:3, 1:3) = R;
    SE(1:3, 4) =  p;
    SE(4, :) = [0 0 0 1];
end

function mult = matrixMulti(m1, m2)
    mult = m1 * m2;
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
