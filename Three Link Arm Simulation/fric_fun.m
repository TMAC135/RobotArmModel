function fric = fric_fun(in1)
%FRIC_FUN
%    FRIC = FRIC_FUN(IN1)

%    This function was generated by the Symbolic Math Toolbox version 5.11.
%    06-Mar-2015 16:41:21

x4 = in1(4,:);
x5 = in1(5,:);
x6 = in1(6,:);
fric = [x4.*4.0;x5.*4.0;x6.*4.0];
