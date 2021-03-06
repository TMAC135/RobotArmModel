function signaling = signaling_fun(in1,in2)
%SIGNALING_FUN
%    SIGNALING = SIGNALING_FUN(IN1,IN2)

%    This function was generated by the Symbolic Math Toolbox version 5.11.
%    06-Mar-2015 16:41:23

l1 = in2(1,:);
l2 = in2(2,:);
l3 = in2(3,:);
x1 = in1(1,:);
x2 = in1(2,:);
x3 = in1(3,:);
x7 = in1(7,:);
x8 = in1(8,:);
x9 = in1(9,:);
t3 = cos(x1);
t4 = cos(x7);
t2 = t3-t4;
t8 = x1+x2;
t9 = x7+x8;
t20 = cos(t8);
t21 = l2.*t20;
t22 = cos(t9);
t23 = l2.*t22;
t24 = l1.*t3;
t25 = l1.*t4;
t5 = t21-t23+t24-t25;
t6 = l1.^2;
t10 = sin(x1);
t11 = sin(x7);
t7 = t10-t11;
t13 = sin(t8);
t14 = l2.*t13;
t15 = sin(t9);
t16 = l2.*t15;
t17 = l1.*t10;
t18 = l1.*t11;
t12 = t14-t16+t17-t18;
t26 = x1+x2+x3;
t27 = x7+x8+x9;
t19 = t14-t16+t17-t18+l3.*sin(t26)-l3.*sin(t27);
t28 = t21-t23+t24-t25+l3.*cos(t26)-l3.*cos(t27);
signaling = t2.^2.*t6+t6.*t7.^2+t5.^2+t12.^2+t19.^2+t28.^2;
