clear all
close all
clc

addpath('E:\College\2021 Fall\CS 403\classCode\matlab_utils');
addpath('E:\College\2021 Fall\CS 403\myCode\hw3');

% each joint consider as a frame
% the position of each frame respect to its last frame
p0 = [0 0 0];
p01 = [0 0 0];
p12 = [0 0 0.15];
p23 = [0.3 0 0];
p34 = [0.15 0 0];
p45 = [0.1 0 0];
p56 = [0.07 0 0];
p6ee = [0.05 0 0]; %end effector

%% Initial Robot figure
% the position of each frame respect to frame 0
f0 = p0;
f1 = f0 + p01;
f2 = f1 + p12;
f3 = f2 + p23;
f4 = f3 + p34;
f5 = f4 + p45;
f6 = f5 + p56;
fee = f6 + p6ee;

figure
grid on
hold on
axis equal

% draw line between each frame
drawLine3D(f0, f1);
drawLine3D(f1, f2);
drawLine3D(f2, f3);
drawLine3D(f3, f4);
drawLine3D(f4, f5);
drawLine3D(f5, f6);
drawLine3D(f6, fee);

% draw coordinate on each frame
% this is the original position
drawCoordinate3DScale(eye(3), f0, 0.05);
drawCoordinate3DScale(eye(3), f1, 0.05);
drawCoordinate3DScale(eye(3), f2, 0.05);
drawCoordinate3DScale(eye(3), f3, 0.05);
drawCoordinate3DScale(eye(3), f4, 0.05);
drawCoordinate3DScale(eye(3), f5, 0.05);
drawCoordinate3DScale(eye(3), f6, 0.05);
drawCoordinate3DScale(eye(3), fee, 0.05);
text(0, 0, 0, 'frame{0}');

hold off

%% 1 a)
th1 = 0;
th2 = pi/2;
th3 = 0;
th4 = pi/6;
th5 = pi/2;
th6 = 0;

%% 1 b)
% th1 = 0;
% th2 = 2*pi/3;
% th3 = 0;
% th4 = pi/3;
% th5 = pi/2;
% th6 = 0;

%% 2 a)
% th1 = 0;
% th2 = pi/2;
% th3 = pi/2;
% th4 = pi/6;
% th5 = pi/2;
% th6 = 0;

%% 2 b)
% th1 = 0;
% th2 = pi/3;
% th3 = pi/4;
% th4 = pi/3;Toee
% th5 = pi/2;
% th6 = 0;

%% Build SE3 for each frame
th = [th1 th2 th3 th4 th5 th6];

R1 = rz(th(1));
R2 = ry(th(2));
R3 = ry(th(3));
R4 = rx(th(4));
R5 = ry(th(5));
R6 = rx(th(6));

% SE3 for each frame wrt its last frame
T01 = SE3(R1, p01);
T12 = SE3(R2, p12);
T23 = SE3(R3, p23);
T34 = SE3(R4, p34);
T45 = SE3(R5, p45);
T56 = SE3(R6, p56);
T6ee = SE3(eye(3), p6ee);

% SE3 wrt frame0
T02 = T01 * T12;
T03 = T02 * T23;
T04 = T03 * T34;
T05 = T04 * T45;
T06 = T05 * T56;
T0ee = T06 * T6ee;


%% Robot Simulation
figure
grid on
hold on
axis equal

drawLine3D(p0, T01(1:3, 4));
drawLine3D(T01(1:3, 4), T02(1:3, 4));
drawLine3D(T02(1:3, 4), T03(1:3, 4));
drawLine3D(T03(1:3, 4), T04(1:3, 4));
drawLine3D(T04(1:3, 4), T05(1:3, 4));
drawLine3D(T05(1:3, 4), T06(1:3, 4));
drawLine3D(T06(1:3, 4), T0ee(1:3, 4));

drawCoordinate3DScale(eye(3), p0, 0.05);
drawCoordinate3DScale(T01(1:3, 1:3), T01(1:3, 4), 0.05);
drawCoordinate3DScale(T02(1:3, 1:3), T02(1:3, 4), 0.05);
drawCoordinate3DScale(T03(1:3, 1:3), T03(1:3, 4), 0.05);
drawCoordinate3DScale(T04(1:3, 1:3), T04(1:3, 4), 0.05);
drawCoordinate3DScale(T05(1:3, 1:3), T05(1:3, 4), 0.05);
drawCoordinate3DScale(T06(1:3, 1:3), T06(1:3, 4), 0.05);
drawCoordinate3DScale(T0ee(1:3, 1:3), T0ee(1:3, 4), 0.05);
text(0, 0, 0, 'frame{0}');

hold off

%% helper functions
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