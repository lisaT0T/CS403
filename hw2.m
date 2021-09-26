%question 4
%x stand for theta 1, and y stand for theta 2
syms x y l_1 l_2 l_3
%the jacobian matrix with respected to x, y 
jacobian([l_1*cos(x)+l_2*cos(x-y)-l_3*cos(x), l_1*sin(x)+l_2*sin(x-y)+l_3*sin(x)], [x, y])
