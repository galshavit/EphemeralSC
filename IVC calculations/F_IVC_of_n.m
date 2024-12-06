%%
clear
close all
clc
tic
%%
n = -1.5 *1e12*linspace(0.97,1.025,50);
% Ds = linspace(23,25,181);
D = 24;  % meV  V/nm --> D*3/200
g = 1.45*1e-11; % meV cm^2
numx=7500;
% numx=2400;
%
Delta = linspace(0,4,61);
temp = zeros(length(n),length(Delta));
mus=temp; Fkin=temp; clear('temp');
%
%%
[E,x] = BandStructure6(numx,D);
ek = squeeze(E(:,:,3)); clear('E');
Omega = 1/(((x(1)-x(2))^2)/(4*pi^2))*1e4;
for w=1:length(n)
    %%
    for i=1:length(Delta)
        [mus(w,i),a,b] = get_mu_Delta(n(w),ek,Delta(i),Omega);
        Fkin(w,i) = (sum(a(a<=mus(w,i))) + sum(b(b<=mus(w,i))))/Omega ;
        disp(i);
        toc
    end
    disp(w);
    toc
end
% save('FIVC')
%%
g = 1.46*1e-11; % meV cm^2
opt_D=0*n;
F = Fkin-Fkin(:,1)+repmat(Delta.^2,length(n),1)/g;
figure
hold on
for w=1:length(n)
    plot(Delta,Fkin(w,:)-Fkin(w,1)+Delta.^2/g,'.-')
    [~,ind]=min(Fkin(w,:)-Fkin(w,1)+Delta.^2/g);
    opt_D(w)=Delta(ind);
end
save('FivcN2');

%%
figure
pcolor(Delta,n*1e-12,F);
shading flat
set(gca, 'Layer','top')
[g,ind,phi_0] = get_g(F,Delta,n);
clim([-1 1]*g)
% xlim(xlims)
% ylim(ylims)
colormap(jet(200))
box on
ax = gca;
ax.XAxis.FontSize = 18;
ax.YAxis.FontSize = 18;
set(gca,'TickLabelInterpreter', 'latex');
xlabel('$\Delta_{\rm IVC}$ [meV]','interpreter','latex','FontSize',18);
ylabel('$n$ [$10^{12}$ cm$^{-2}$]','interpreter','latex','FontSize',18);
col=colorbar;
col.FontSize=16;
col.TickLabelInterpreter= 'latex';