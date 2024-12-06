%%
clear
close all
clc
load('FIVC.mat','Ds','Delta','Fkin','n','g')
Delta_test = 2.935;
Ds_jump = 24.145;
%%
numx=2800;
eps=linspace(-50,-25,901);
Ds_demo = linspace(Ds_jump-3,Ds_jump+1,81);
omegad=1;
gatt = 1.18e-12;
%%
de=eps(2)-eps(1);
dos_demo=zeros(length(eps),length(Ds_demo));
dos_demo0=zeros(length(eps),length(Ds_demo));
mu = 0*Ds_demo; mu0 = 0*Ds_demo;
Tc = 0*Ds_demo; Tc0 = 0*Ds_demo;
%%
for w=1:length(Ds_demo)
    %%
    [E,x] = BandStructure6(numx,Ds_demo(w));
    ek = squeeze(E(:,:,3)); clear('E');
    Omega = 1/(((x(1)-x(2))^2)/(4*pi^2))*1e4;
    [mu(w),a,b] = get_mu_Delta(n,ek,Delta_test,Omega);
    [mu0(w),~,~] = get_mu_Delta(n,ek,0,Omega);
    for i=2:length(eps)
        dos_demo(i,w) = sum (a>=eps(i-1) & a<eps(i),'all') + sum (b>=eps(i-1) & b<eps(i),'all');
        dos_demo0(i,w) = 2*sum (ek>=eps(i-1) & ek<eps(i),'all');
    end
    dos_demo(:,w) = dos_demo(:,w)/de/Omega; dos_demo0(:,w) = dos_demo0(:,w)/de/Omega;
    %%
    [Tc(w),~] = get_Tc(eps,dos_demo(:,w),mu(w),omegad,gatt,1e-2);
    [Tc0(w),~] = get_Tc(eps,dos_demo0(:,w),mu(w),omegad,gatt,1e-2);
    %%
    disp(w);
end
%%
figure; hold on;
Tdisplay = Tc0; Tdisplay(Ds_demo>=Ds_jump)=Tc(Ds_demo>=Ds_jump);
plot(Ds_demo*3/200,Tc0*11.6*1e3,'ok','MarkerSize',4)
plot(Ds_demo*3/200,Tdisplay*11.6*1e3,'.','Color',[0.2 0.2 1],'MarkerSize',14)
plot([1 1]*Ds_jump*3/200*0.999,[0 120],'-','Color',[1 1 1]*0.5,'LineWidth',2)
box on
ax = gca;
ax.XAxis.FontSize = 18;
ax.YAxis.FontSize = 18;
set(gca,'TickLabelInterpreter', 'latex');
xlabel('$D$ [V/nm]','interpreter','latex','FontSize',18);
ylabel('$T_c$ [mK]','interpreter','latex','FontSize',18);
ylim([0 110])
xlim([-inf 0.37])
%%
w=30;
figure
hold on
plot(eps-mu0(w),dos_demo0(:,w),'.-')
plot(eps-mu(w),dos_demo(:,w),'.-')
ylim([0 inf])
xlim([-10 8])
plot([1 1]*omegad, 8.2e11*[0 1],'--k')
plot([-1 -1]*omegad, 8.2e11*[0 1],'--k')
box on
ax = gca;
ax.XAxis.FontSize = 18;
ax.YAxis.FontSize = 18;
set(gca,'TickLabelInterpreter', 'latex');
xlabel('$\xi$ [meV]','interpreter','latex','FontSize',18);
ylabel('${\cal N}$ [meV$^{-1}$cm$^{-2}$]','interpreter','latex','FontSize',18);
leg=legend('Normal','IVC','','');
leg.FontSize=16;
leg.Interpreter= 'latex';