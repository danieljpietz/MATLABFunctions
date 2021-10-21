function [H,d, G] = linkDynamics(T,JJ,GAMMA,J,dotJ,m, dotgamma, gg)
temp=J*dotgamma;
w=temp(1:3);
M=[JJ,skew(GAMMA)*T.';T*skew(GAMMA).',m*eye(3)];
H=J.'*M*J;
d=J.'*M*dotJ*dotgamma + J.'*[cross(w,JJ*w);T*cross(w,cross(w,GAMMA))];
G=-J.'*[skew(GAMMA)*T.'*gg;m*gg];
end

