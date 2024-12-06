%%
% in units where
% 16g = 1
% \xi_\phi = 1
% thus, \sigma=0.5
clear
close all
clc
%%
xi_psi = 20;
a0 = 0.005;
alpha = 0.006;
B=linspace(0,0.1,501);
Jdw=0*B; Rc=0*B; fSC=0*B; R0=0*B;
%%
for j=1:length(B)
    Jdw(j) = DW_energy(B(j));
    fSC(j) = (a0+alpha*B(j))^2/(4*a0);
    Rc(j) = (Jdw(j)+8/3*xi_psi*fSC(j)+fSC(j)*0.5*(1+xi_psi)) / (B(j)-fSC(j));
    R0(j) = (Jdw(j)) / (B(j));
end
E = pi*Rc.^2.*(-B) + pi*(Rc+0.5*(1+xi_psi)).^2.*fSC + 2*pi*Rc.*Jdw + 2*pi*(Rc+0.5*(1+xi_psi)).*8/3.*xi_psi.*fSC;
E0 = pi*R0.^2.*(-B)  + 2*pi*R0.*Jdw;
Rc(fSC>B)=nan;
E(fSC>B)=nan;
%%
figure
hold on
plot(B,Rc,'.-')
hold on
plot(B,R0,'--','Color',[0.5 0.1 0.1],'LineWidth',2)
set(gca,'yscale','log')
box on
ax = gca;
ax.XAxis.FontSize = 18;
ax.YAxis.FontSize = 18;
set(gca,'TickLabelInterpreter', 'latex');
xlabel('$B/\left(16g\right)$','interpreter','latex','FontSize',18);
ylabel('$R_c/\xi_\phi$','interpreter','latex','FontSize',18);
ylim([1 inf])
yticks([1,10,100,1000])
%%
figure
hold on
plot(B,E,'.-')
hold on
plot(B,E0,'--','Color',[0.5 0.1 0.1],'LineWidth',2)
set(gca,'yscale','log')
box on
ax = gca;
ax.XAxis.FontSize = 18;
ax.YAxis.FontSize = 18;
set(gca,'TickLabelInterpreter', 'latex');
xlabel('$B/\left(16g\right)$','interpreter','latex','FontSize',18);
ylabel('$E_{\rm barrier}/\left(16g\xi_\phi^2\right)$','interpreter','latex','FontSize',18);
ylim([1 inf])
