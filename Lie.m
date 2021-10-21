function [LpQ] = Lie(P, Q, degree)
    if nargin ~= 3
        degree = 1;
    end
    if degree == 1
        LpQ = sum(gradient(Q).*P);
    else
        LpQ = Lie(P, Lie(P, Q, degree - 1));
    end
end

