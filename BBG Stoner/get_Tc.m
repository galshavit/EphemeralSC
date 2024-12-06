function [Tc,Fcond] = get_Tc(eps,dos,EF,omegad,gatt,Tprobe)
Ts = linspace(1,500,500)/(11.6*1e3); calD=0*Ts;
Delta=1.76*Ts; F=0*Delta;

epss = linspace(min(eps),max(eps),length(eps)*10);
dos = interp1(eps,dos,epss);
eps=epss;

de = eps(2)-eps(1);
xi = eps-EF;
cond = abs(xi)<omegad;

U0 = 1 / (2*interp1(eps,dos,EF));
ell = 2*de*sum(~cond.*dos./abs(xi));
geff = gatt - U0/(1+U0*ell);

for j=1:length(Ts);     calD(j) = 2*de*sum(cond.*dos.*tanh(xi/(2*Ts(j)))./xi);    end
Tc = interp1(geff*calD+cumsum(calD*0+1e-8),Ts,1);

for j=1:length(Delta)
    F(j) = 0.5*Delta(j)^2/geff...
        - 2*de*sum(cond.*dos.*(sqrt(xi.^2+Delta(j)^2)-abs(xi)))...
        -4*Tprobe*de*sum(cond.*dos.*log((1+exp(-sqrt(xi.^2+Delta(j)^2)/Tprobe))./(1+exp(-abs(xi)/Tprobe))));
end
Fcond = min(F);
end