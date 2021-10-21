function V = sym2poly2(P, var)
syms f(t)
V = [];


while f(var) ~= 0
   V = [V, f(0)];
   f(var) = diff(f, var);
end

end