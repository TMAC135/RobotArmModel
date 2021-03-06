function C_jac = C_jac_fun(in1)
%C_JAC_FUN
%    C_JAC = C_JAC_FUN(IN1)

%    This function was generated by the Symbolic Math Toolbox version 5.11.
%    06-Mar-2015 16:41:22

x2 = in1(2,:);
x3 = in1(3,:);
x4 = in1(4,:);
x5 = in1(5,:);
x6 = in1(6,:);
t2 = x2+x3;
t3 = cos(t2);
t4 = t3.*(9.0./2.0e2);
t5 = sin(t2);
t6 = cos(x2);
t7 = t6.*(1.4e1./2.5e1);
t8 = t4+t7;
t9 = cos(x3);
t10 = t9.*(9.0./2.5e2);
t11 = t4+t10;
t12 = sin(x2);
t14 = t5.*(9.0./2.0e2);
t15 = t12.*(1.4e1./2.5e1);
t13 = -t14-t15;
t16 = sin(x3);
t17 = x4+x5+x6;
t19 = t16.*(9.0./2.5e2);
t18 = -t14-t19;
t20 = t3.*x4.*(9.0./2.0e2);
C_jac = reshape([0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t3.*x6.*(-9.0./2.0e2)-t8.*x5,t8.*x4,t20,t3.*x6.*(-9.0./2.0e2)-t8.*x4-t8.*x5,0.0,0.0,t3.*t17.*(-9.0./2.0e2),0.0,0.0,t3.*x5.*(-9.0./2.0e2)-t11.*x6,t20-t9.*x6.*(9.0./2.5e2),t9.*x5.*(9.0./2.5e2)+t11.*x4,t3.*x4.*(-9.0./2.0e2)-t3.*x5.*(9.0./2.0e2)-t11.*x6,t9.*x6.*(-9.0./2.5e2),t9.*(x4+x5).*(9.0./2.5e2),t17.*(t3.*5.0+t9.*4.0).*(-9.0./1.0e3),t9.*t17.*(-9.0./2.5e2),0.0,0.0,t14+t15,t14+t19,t13,0.0,t19,t18,-t19,0.0,t13,0.0,t19,t13,0.0,t19,t18,-t19,0.0,t5.*(-9.0./2.0e2)-t16.*(9.0./2.5e2),-t19,0.0,t18,-t19,0.0,t18,-t19,0.0],[3, 3, 6]);
