function [Jdw] = DW_energy(B)
% in units where
% 16g = 1
% \xi_\phi = 1
% thus, \sigma=0.5
phi = linspace(0,1-sqrt(B)-1e-16,1e4+1);
dphi = phi(2)-phi(1);
Jdw = dphi*sum(sqrt(phi.^2.*(phi-1).^2-B*phi.^2));
end