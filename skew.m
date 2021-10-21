function [S] = skew(wOld)
S = [0 -wOld(3) wOld(2) ; wOld(3) 0 -wOld(1) ; -wOld(2) wOld(1) 0 ];
end

