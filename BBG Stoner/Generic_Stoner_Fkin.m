function [Fkin] = Generic_Stoner_Fkin(eps,dos,n_eps,ntot,d,delta)
%
E0 = interp1(n_eps,eps,ntot/(2*d));
Ep = interp1(n_eps,eps,ntot/(2*d)+delta);
Em = interp1(n_eps,eps,ntot/(2*d)-delta);
%
cond1 = eps<=Ep & eps>=E0;
cond2 = eps<=E0 & eps>=Em;
de = eps(2)-eps(1);
%
Fkin = d*de*sum((eps-E0).*dos.*cond1 - (eps-E0).*dos.*cond2);
end