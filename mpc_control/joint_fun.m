function joint = joint_fun(in1)
%JOINT_FUN
%    JOINT = JOINT_FUN(IN1)

%    This function was generated by the Symbolic Math Toolbox version 5.10.
%    07-Jan-2016 21:13:43

x1 = in1(1,:);
joint = [cos(x1);sin(x1)];
