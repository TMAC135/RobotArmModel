function signaling_jac2 = signaling_jac2_fun(in1,in2)
%SIGNALING_JAC2_FUN
%    SIGNALING_JAC2 = SIGNALING_JAC2_FUN(IN1,IN2)

%    This function was generated by the Symbolic Math Toolbox version 5.11.
%    06-Mar-2015 16:41:26

l1 = in2(1,:);
l2 = in2(2,:);
l3 = in2(3,:);
x1 = in1(1,:);
x2 = in1(2,:);
x3 = in1(3,:);
x7 = in1(7,:);
x8 = in1(8,:);
x9 = in1(9,:);
t2 = l2.^2;
t3 = x1+x2-x7-x8;
t4 = cos(t3);
t5 = t2.*t4.*4.0;
t6 = l3.^2;
t7 = x1+x2+x3-x7-x8-x9;
t8 = cos(t7);
t9 = t6.*t8.*2.0;
t10 = x1+x2+x3-x7-x8;
t11 = cos(t10);
t12 = l2.*l3.*t11.*2.0;
t13 = x1+x2+x3-x7;
t14 = cos(t13);
t15 = l1.*l3.*t14.*2.0;
t16 = x1+x2-x7-x8-x9;
t17 = cos(t16);
t18 = l2.*l3.*t17.*2.0;
t19 = x1+x2-x7;
t20 = cos(t19);
t21 = l1.*l2.*t20.*4.0;
t22 = l1.^2;
t23 = x1-x7;
t24 = cos(t23);
t25 = t22.*t24.*6.0;
t26 = x1-x7-x8;
t27 = cos(t26);
t28 = l1.*l2.*t27.*4.0;
t29 = x1-x7-x8-x9;
t30 = cos(t29);
t31 = l1.*l3.*t30.*2.0;
t32 = t5+t9+t12+t15+t18+t21;
t33 = x2+x3;
t34 = cos(t33);
t35 = t9+t12+t15;
t37 = l1.*l3.*t34.*2.0;
t36 = t9+t12+t15-t37;
t38 = -t5-t9-t12-t15-t18-t21-t25-t28-t31;
t39 = -t5-t9-t12-t15-t18-t21;
t40 = -t9-t12-t15;
t41 = t5+t9+t12+t15+t18+t21+t25+t28+t31;
t42 = -t5-t9-t12-t18-t28-t31;
t43 = -t5-t9-t12-t18;
t44 = -t9-t12;
t45 = t5+t9+t12+t18+t28+t31;
t46 = x8+x9;
t47 = cos(t46);
t48 = -t9-t18-t31;
t49 = -t9-t18;
t50 = t9+t18+t31;
t52 = l1.*l3.*t47.*2.0;
t51 = t9+t18+t31-t52;
signaling_jac2 = reshape([t41,t32,t35,0.0,0.0,0.0,t38,t42,t48,0.0,0.0,0.0,t32,t5+t9+t12+t15+t18+t21-l1.*l2.*cos(x2).*4.0-l1.*l3.*t34.*2.0,t36,0.0,0.0,0.0,t39,t43,t49,0.0,0.0,0.0,t35,t36,t9+t12+t15-t37-l2.*l3.*cos(x3).*2.0,0.0,0.0,0.0,t40,t44,-t9,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t38,t39,t40,0.0,0.0,0.0,t41,t45,t50,0.0,0.0,0.0,t42,t43,t44,0.0,0.0,0.0,t45,t5+t9+t12+t18+t28+t31-l1.*l2.*cos(x8).*4.0-l1.*l3.*t47.*2.0,t51,0.0,0.0,0.0,t48,t49,-t9,0.0,0.0,0.0,t50,t51,t9+t18+t31-t52-l2.*l3.*cos(x9).*2.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0],[12, 12]);
