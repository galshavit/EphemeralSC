%%
clear
close all
clc
load('RTG_DOS.mat'); Ds_dos=Ds;
load('RTG_IVC_data.mat')
%%
dos_ivc = zeros(length(Ds),length(eps));
for i=1:length(eps)
    dos_ivc(:,i) = interp1(Ds_dos,dos(:,i),Ds);
end
dos=dos_ivc;
%%
positive_threshold = 0;
xi_phi = 1/sqrt(abs(ntot));
gatt=2.1e-12;
omegad=0.5;
clims = [-1 1]*2e9;
xlims = [0 3.8];
ylims = [0.346 0.374];
Jdw = 0*Ds;
B=0*Ds;
EF = 0*Ds;
Tc=0*Ds;
Fcond=0*Ds;
%% Plot Free Energy Landscape 
figure
pcolor(Delta,Ds*3/200,F);
shading flat
set(gca, 'Layer','top')
[g,ind,phi_0] = get_g(F,Delta,Ds);
clim([-1 1]*g)
xlim(xlims)
ylim(ylims)
colormap(jet(200))
box on
ax = gca;
ax.XAxis.FontSize = 18;
ax.YAxis.FontSize = 18;
set(gca,'TickLabelInterpreter', 'latex');
xlabel('$\Delta_{\rm IVC}$ [meV]','interpreter','latex','FontSize',18);
ylabel('D [V/nm]','interpreter','latex','FontSize',18);
col=colorbar;
col.FontSize=16;
col.TickLabelInterpreter= 'latex';
%%
figure
plot(Delta,F(74,:),'o-','Color',[0.6 0.11 0.5],'LineWidth',1.5,'MarkerSize',5)
box on
ax = gca;
ax.XAxis.FontSize = 18;
ax.YAxis.FontSize = 18;
set(gca,'TickLabelInterpreter', 'latex');
xlabel('$\Delta_{\rm IVC}$ [meV]','interpreter','latex','FontSize',18);
ylabel('$F$ [meV/cm$^2$]','interpreter','latex','FontSize',18);
ylim([-0.7,4]*1e9)
xlim(xlims)
%% Calculate independent Energy barrier and SC Tc
for i=1:length(Ds)
    f = F(i,:);
    Jdw(i) = J_phi(Delta, f , g, xi_phi, positive_threshold,phi_0);
    B(i) = -min(f);
    EF(i) = interp1(n_of_eps(eps,dos(i,:))+cumsum(eps*0+1),eps,ntot/(2*2));
    [Tc(i),Fcond(i)] = get_Tc(eps,dos(i,:),EF(i),omegad,gatt,20/(11.6*1e3));
end
figure
hold on
plot(Ds*3/200,pi*(Jdw.^2)./B,'.')
% plot(276/(11.6*1e3)+Ds*0,Ds*3/200,'--k')
box on
ax = gca;
ax.XAxis.FontSize = 18;
ax.YAxis.FontSize = 18;
set(gca,'TickLabelInterpreter', 'latex');
ylabel('$E_{\rm barrier}$ [meV]','interpreter','latex','FontSize',18);
xlabel('D [V/nm]','interpreter','latex','FontSize',18);
xlim(ylims)
% set(gca,'xscale','log')
%% Time for superconductivity
figure
plot(Tc*11.6*1e3,Ds*3/200,'.-')
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
plot(Ds*3/200,E*11.6*1e3/Tmk,'.-')
hold on
plot(Ds*3/200,pi*(Jdw.^2)./B*11.6*1e3/Tmk,'--','Color',[0.5 0.1 0.1],'LineWidth',2)
set(gca,'Yscale', 'log');
xlim([0.356 0.366])
ylim([1 500])
box on
ax = gca;
ax.XAxis.FontSize = 18;
ax.YAxis.FontSize = 18;
set(gca,'TickLabelInterpreter', 'latex');
ylabel('$\beta E_{\rm barrier}$','interpreter','latex','FontSize',18);
xlabel('D [V/nm]','interpreter','latex','FontSize',18);

figure
hold on
plot(Ds*3/200,Jpsi./Jdw,'.-','Color',[0.2 0.6 0.1])
ylim([0 1])
xlim([0.356 0.366])

box on
ax = gca;
ax.XAxis.FontSize = 18;
ax.YAxis.FontSize = 18;
set(gca,'TickLabelInterpreter', 'latex');
ylabel('$J_\phi/J_\Psi$ ','interpreter','latex','FontSize',18);
xlabel('D [V/nm]','interpreter','latex','FontSize',18);
xlim([0.356 0.366])
%%
figure
hold on
ax = nexttile;
yyaxis(ax,'left')
plot(ax,Ds*3/200,Tc*11.6*1e3,'.-')
ylabel('$T_c$ [mK] ','interpreter','latex','FontSize',18);
yyaxis(ax,'right')
plot(Ds*3/200,Fcond,'.-')
ylabel('$F_{\rm cond}$ [meV cm$^{-2}$] ','interpreter','latex','FontSize',18);
box on
ax = gca;
ax.XAxis.FontSize = 18;
set(ax(1), 'FontSize', 18)
set(gca,'TickLabelInterpreter', 'latex');
xlabel('D [V/nm]','interpreter','latex','FontSize',18);
xlim([-inf inf])
%%
figure
plot(Ds*3/200,2*Rc*1e6,'.-','Color',[0.1 0 0.5])
% set(gca,'Yscale', 'log');
box on
ax = gca;
ax.XAxis.FontSize = 18;
ax.YAxis.FontSize = 18;
set(gca,'TickLabelInterpreter', 'latex');
ylabel('$2R_c$ [$\mu$m]','interpreter','latex','FontSize',18);
xlabel('D [V/nm]','interpreter','latex','FontSize',18);
xlim([0.356 0.366])
ylim([0 10])