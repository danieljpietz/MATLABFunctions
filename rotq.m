function out1 = rotq(q)
%ROTQOPTIMIZED
%    OUT1 = ROTQOPTIMIZED(q(1),q(2),q(3),q(4))

%    This function was generated by the Symbolic Math Toolbox qersion 8.4.
%    01-May-2020 16:12:25


t2 = q(1).^2;
t3 = q(2).^2;
t4 = q(3).^2;
t5 = q(4).^2;
t6 = q(1).*q(2).*2.0;
t7 = q(1).*q(3).*2.0;
t8 = q(1).*q(4).*2.0;
t9 = q(2).*q(3).*2.0;
t10 = q(2).*q(4).*2.0;
t11 = q(3).*q(4).*2.0;
t12 = -t3;
t13 = -t4;
t14 = -t5;
out1 = reshape([t2+t3+t13+t14,-t8+t9,t7+t10,t8+t9,t2+t4+t12+t14,-t6+t11,-t7+t10,t6+t11,t2+t5+t12+t13],[3,3]);
