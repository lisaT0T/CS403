function gravity = grav_pend(in1,in2)
%GRAV_PEND
%    GRAVITY = GRAV_PEND(IN1,IN2)

%    This function was generated by the Symbolic Math Toolbox version 8.7.
%    12-Nov-2021 15:07:43

c1 = in2(1,:);
c2 = in2(3,:);
g = in2(9,:);
l1 = in2(2,:);
m1 = in2(5,:);
m2 = in2(7,:);
th1 = in1(1,:);
th2 = in1(2,:);
t2 = sin(th1);
t3 = th1+th2;
t4 = sin(t3);
gravity = [g.*m2.*(c2.*t4+l1.*t2)+c1.*g.*m1.*t2;c2.*g.*m2.*t4];