function phi = phi_fun(in1)
%PHI_FUN
%    PHI = PHI_FUN(IN1)

%    This function was generated by the Symbolic Math Toolbox version 5.10.
%    07-Jan-2016 20:46:56

k1 = in1(3,:);
k2 = in1(4,:);
q = in1(1,:);
t2 = sin(q);
t3 = cos(q);
phi = k1.*t2.*(2.0./5.0)-k2.*t2.*(2.0./5.0)-k1.*t2.*1.0./sqrt(t3.*-2.0e1+2.9e1).*(8.0./5.0)+k2.*t2.*1.0./sqrt(t3.*2.0e1+2.9e1).*(8.0./5.0);
