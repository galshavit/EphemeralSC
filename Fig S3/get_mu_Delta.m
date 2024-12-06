function [EF,a,b] = get_mu_Delta(ntot,ek,Delta,Omega)
eps=linspace(max(ek,[],'all')-10,max(ek,[],'all'),1501); ns=0*eps;
%%
ek2=flipud(ek);
E1 = 0.5*(ek+ek2) + sqrt(Delta^2+0.25*(ek-ek2).^2);
E2 = 0.5*(ek+ek2) - sqrt(Delta^2+0.25*(ek-ek2).^2);
a=sort(E1(:)); b=sort(E2(:));
%%
for i=1:length(eps)
    ns(i) = -sum((a>=eps(i)),'all')/Omega -sum((b>=eps(i)),'all')/Omega;
end

EF = interp1(ns+cumsum(0*ns+1),eps,ntot/2);
end