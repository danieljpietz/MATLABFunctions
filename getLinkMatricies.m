function [H_LINK,d_LINK, G_LINK, J_LINK, dotJ_LINK] = getLinkMatricies(T_LINK,r_II_LINK, JJ_LINK, m_LINK, GAMMA_LINK, gamma, dotgamma, g_dir,g)
if ~exist('g','var')
   g = 9.81;
end
J_LINK=simplify([T_LINK(:,3).'*jacobian(T_LINK(:,2),gamma);...
    T_LINK(:,1).'*jacobian(T_LINK(:,3),gamma);...
    T_LINK(:,2).'*jacobian(T_LINK(:,1),gamma);...
    jacobian(r_II_LINK,gamma)]);
dotJ_LINK=sym(zeros(6, length(gamma)));
for i = 1:length(gamma)
    dotJ_LINK(:,i) = jacobian(J_LINK(:,i),gamma)*dotgamma;
end
dotJ_LINK = simplify(dotJ_LINK);
w_LINK=J_LINK(1:3,:)*dotgamma;

M_LINK=[JJ_LINK,skew(GAMMA_LINK)*T_LINK.';...
    T_LINK*skew(GAMMA_LINK).',m_LINK*eye(3)];
H_LINK=simplify(J_LINK.'*M_LINK*J_LINK);
d_LINK=simplify(J_LINK.'*M_LINK*dotJ_LINK*dotgamma + J_LINK.'*[cross(w_LINK,JJ_LINK*w_LINK);T_LINK*cross(w_LINK,cross(w_LINK,GAMMA_LINK))]);
G_LINK=simplify(-J_LINK.'*[skew(GAMMA_LINK)*T_LINK.'*g_dir;m_LINK*g_dir]);

end

