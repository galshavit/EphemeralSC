%%
clear
close all
clc
load('BBG_Stoner_data.mat')
%%
positive_threshold = 0;
xi_phi = 1/sqrt(abs(ntot));
gatt=2.65e-12;
omegad=0.5;
clims = [-1 1]*2e9;
xlims = [0 13.5]*1e10;
ylims = [0.98 1.03];
Jdw = 0*Ds;
B=0*Ds;
EF = 0*Ds;
Tc=0*Ds;
Fcond=0*Ds;
%% Plot Free Energy Landscape 
figure
pcolor(delta,Ds,F);
shading interp
set(gca, 'Layer','top')
[g,ind,phi_0] = get_g(F,delta,Ds);
clim([-1 1]*g)
xlim(xlims)
ylim(ylims)
colormap(jet(200))
box on
ax = gca;
ax.XAxis.FontSize = 18;
ax.YAxis.FontSize = 18;
set(gca,'TickLabelInterpreter', 'latex');
xlabel('$\delta$ [cm$^{-2}$]','interpreter','latex','FontSize',18);
ylabel('D [V/nm]','interpreter','latex','FontSize',18);
col=colorbar;
col.FontSize=16;
col.TickLabelInterpreter= 'latex';
figure
hold on
plot(delta,F(325,:),'o-','Color',[0.6 0.11 0.5],'LineWidth',0.5,'MarkerSize',5)
box on
ax = gca;
ax.XAxis.FontSize = 18;
ax.YAxis.FontSize = 18;
set(gca,'TickLabelInterpreter', 'latex');
xlabel('$\delta$ [cm$^{-2}$]','interpreter','latex','FontSize',18);
ylabel('$F$ [meV/cm$^2$]','interpreter','latex','FontSize',18);
ylim([-0.95,4]*1e9)
xlim([0 1.4e11])
%% Calculate independent Energy barrier and SC Tc
for i=1:length(Ds)
    f = F(i,:);
    Jdw(i) = J_phi(delta, f , g, xi_phi, positive_threshold,phi_0);
    B(i) = -min(f);
    EF(i) = interp1(n_of_eps(eps,dos(i,:))+cumsum(eps*0+1),eps,ntot/(2*2));
    [Tc(i),Fcond(i)] = get_Tc(eps,dos(i,:),EF(i),omegad,gatt,30/(11.6*1e3));
end
figure
hold on
plot(pi*(Jdw.^2)./B,Ds,'.-')
% plot(pi*(Jdw.^3)./B.^2*20/xi_phi,Ds*3/200,'.-')
plot(200/(11.6*1e3)+Ds*0,Ds,'--k')
ylim(ylims)
xlim([1/(11.6*1e3) inf])
box on
ax = gca;
ax.XAxis.FontSize = 18;
ax.YAxis.FontSize = 18;
set(gca,'TickLabelInterpreter', 'latex');
xlabel('$E_{\rm barrier}$ [meV]','interpreter','latex','FontSize',18);
ylabel('D [V/nm]','interpreter','latex','FontSize',18);
% set(gca,'xscale','log')
%% Time for superconductivity
figure
plot(Tc*11.6*1e3,Ds,'.-')
ylim(ylims)
% xlim([1/(11.6*1e3) inf])
box on
ax = gca;
ax.XAxis.FontSize = 18;
ax.YAxis.FontSize = 18;
set(gca,'TickLabelInterpreter', 'latex');
xlabel('$T_c$ [mK]','interpreter','latex','FontSize',18);
ylabel('D [V/nm]','interpreter','latex','FontSize',18);
%% Proper instanton
Bthresh=4e7;
B(B<Bthresh)=nan;
xi_0 = 150*1e-7;
Tc_0 = 100;
Tmk=30;
xi_psi = xi_0*sqrt(Tc_0^2./(Tc*11.6*1e3)./(Tc*11.6*1e3-Tmk));
Jpsi = 8/3*abs(Fcond).*xi_psi;
Rc = (Jdw+Jpsi+0.5*(xi_phi+xi_psi).*abs(Fcond))./(B+abs(Fcond));
E = -pi*Rc.^2.*B + pi*(Rc+0.5*(xi_phi+xi_psi)).^2.*abs(Fcond)...
    +2*pi*Rc.*Jdw + 2*pi*(Rc+0.5*(xi_phi+xi_psi)).*Jpsi;
E(Tc*11.6*1e3<Tmk)=nan;
plot(Ds,E*11.6*1e3/Tmk,'.-')
hold on
plot(Ds,pi*(Jdw.^2)./B*11.6*1e3/Tmk,'--','Color',[0.5 0.1 0.1],'LineWidth',2)
set(gca,'Yscale', 'log');
xlim([0.998 1.013])
ylim([1 300])
box on
ax = gca;
ax.XAxis.FontSize = 18;
ax.YAxis.FontSize = 18;
set(gca,'TickLabelInterpreter', 'latex');
ylabel('$\beta E_{\rm barrier}$','interpreter','latex','FontSize',18);
xlabel('D [V/nm]','interpreter','latex','FontSize',18);

Jpsi(Tc*11.6*1e3<Tmk)=nan;
Jdw(Jdw<10)=nan;
figure
hold on
plot(Ds,Jpsi./Jdw,'.-','Color',[0.2 0.6 0.1])
ylim([0 0.125])
xlim([0.998 1.013])

box on
ax = gca;
ax.XAxis.FontSize = 18;
ax.YAxis.FontSize = 18;
set(gca,'TickLabelInterpreter', 'latex');
ylabel('$J_\Psi/J_\phi$ ','interpreter','latex','FontSize',18);
xlabel('D [V/nm]','interpreter','latex','FontSize',18);
%%
figure
hold on
ax = nexttile;
yyaxis(ax,'left')
plot(ax,Ds,Tc*11.6*1e3,'.-')
ylabel('$T_c$ [mK] ','interpreter','latex','FontSize',18);
yyaxis(ax,'right')
plot(Ds,Fcond,'.-')
ylabel('$F_{\rm cond}$ [meV cm$^{-2}$] ','interpreter','latex','FontSize',18);
box on
ax = gca;
ax.XAxis.FontSize = 18;
set(ax(1), 'FontSize', 18)
set(gca,'TickLabelInterpreter', 'latex');
xlabel('D [V/nm]','interpreter','latex','FontSize',18);
xlim([0.96 inf])
%%
gamma= (xi_phi/100/(1e5) * exp(15))^(-1);
gamma0= (xi_phi/100/(1e5) * exp(8))^(-1);
gamma_psi = (100*1e-12)^-1;
t = logspace(-20,1,211);
p = gamma_psi/(gamma_psi+gamma0-gamma)*...
    (1-exp(-(gamma_psi+gamma0-gamma)*t)).*exp(-gamma*t);
figure
hold on
plot(t,p,'.-')
% plot(t,exp(-gamma*t))
plot([1 1]/gamma_psi, [0,1],'--')
plot([1 1]/gamma, [0,1],'--')
xlim([1e-13,1e-5])
ylim([0 1])
set(gca,'Xscale', 'log');
box on
ax = gca;
ax.XAxis.FontSize = 12;
ax.YAxis.FontSize = 18;
set(gca,'TickLabelInterpreter', 'latex');
ylabel('$p_{\rm false}$ ','interpreter','latex','FontSize',18);
xlabel('$t$ [sec]' ,'interpreter','latex','FontSize',18);
%%
figure
plot(Ds,2*real(Rc)*1e6,'.-','Color',[0.1 0 0.5])
% set(gca,'Yscale', 'log');
box on
ax = gca;
ax.XAxis.FontSize = 18;
ax.YAxis.FontSize = 18;
set(gca,'TickLabelInterpreter', 'latex');
ylabel('$2R_c$ [$\mu$m]','interpreter','latex','FontSize',18);
xlabel('D [V/nm]','interpreter','latex','FontSize',18);
xlim([0.998 1.013])
ylim([0 10])