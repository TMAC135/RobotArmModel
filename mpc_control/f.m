function fx=f(x,u)
%fx=[x(3); x(4); qddot_fun(x,u)]; %This funciton is for 2_link_arm, needed to be enabled running demo.m
fx=[x(2);qddot_fun(x,u)]; %This function is for the pendulum function, needed
%to be enabled running demo_pendulum.m