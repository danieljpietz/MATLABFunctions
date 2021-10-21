function [TN, wN, JN, dotJN] = recursive_jacobian(T_prev, J_prev, dotJ_prev, T_local, r_local, I_hat, I_tilde, gamma, dotgamma, options)
w_prev = J_prev(1:3,:) * dotgamma;
TN = T_prev * T_local;
w_local = I_hat * dotgamma;
dotr_local = I_tilde * dotgamma;

JN = [T_local.', zeros(3,3);
    -T_prev*skew(r_local), eye(3)] * J_prev + [I_hat;T_prev*I_tilde];
wN = JN(1:3,:)*dotgamma;

dotJN = [-skew(w_local)*T_local.', zeros(3,3); ...
    -T_prev*(skew(w_prev)*skew(r_local)+skew(dotr_local)), zeros(3,3)]*J_prev + ...
    [T_local.', zeros(3,3);
    -T_prev*skew(r_local), eye(3)] * dotJ_prev + [zeros(3,length(gamma));T_prev*skew(w_prev)*I_tilde];

if exist('options','var')
    if strcmp(options, 'Simple')
        TN = simplify(TN);
        wN = simplify(wN);
        JN = simplify(JN);
        dotJN = simplify(dotJN);
    end
    
end


    