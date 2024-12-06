function [Jdw] = J_phi(phi, F_of_phi , g, xi_phi, positive_threshold,phi_0)
dphi = phi(2)-phi(1);
[~,minind] = min(F_of_phi);
cond = ( phi<=phi(minind) ) & ( F_of_phi>=positive_threshold ) ;
Jdw = 4*g*xi_phi *dphi*sum(cond.* sqrt(F_of_phi))/sqrt(g)/phi_0;
end