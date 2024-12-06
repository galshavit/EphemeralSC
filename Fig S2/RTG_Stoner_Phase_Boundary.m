%%
clear
close all
clc
load('RTG_DOS_basin.mat')
tic
%%
n_list = linspace(-12*1e11,-4*1e11,61);
% n_list=-6.2*1e11;
U = 2.1e-11;
d=2;
delta = linspace(0,abs(12*1e11)/(2*d),800);
%%
temp = zeros(length(n_list),length(Ds),length(delta));
Fkin = temp;
F = temp; 
Delta = repmat(delta,length(Ds),1);
%%
for z=1:length(n_list)
    for i=1:length(Ds)
        n_eps = n_of_eps(eps,dos(i,:));
        for j=1:length(delta)
            Fkin(z,i,j)=Generic_Stoner_Fkin(eps,dos(i,:),n_eps + cumsum(n_eps*0+0.01),n_list(z),d,delta(j));
        end
        disp(i);
        toc
    end
    F(z,:,:) = squeeze(Fkin(z,:,:)) - U*Delta.^2;
end
%%
opt_D = zeros(length(n_list),length(Ds));
for z=1:length(n_list)
    for w=1:length(Ds)
        tempo = squeeze(F(z,w,:))'.*(delta<=abs(n_list(z))/(2*d));
        [~,ind]=min(tempo);
        opt_D(z,w)=delta(ind);
    end
end
%%
figure
imagesc(n_list,Ds*3/200,opt_D');
set(gca,'ydir','normal');
ylim([-inf inf])
xlim([-1.1e12 inf])
clim([0.2 1]*1e11)
box on
ax = gca;
ax.XAxis.FontSize = 18;
ax.YAxis.FontSize = 18;
set(gca,'TickLabelInterpreter', 'latex');
xlabel('$n_{\rm tot}$ [cm$^{-2}$]','interpreter','latex','FontSize',18);
ylabel('D [V/nm]','interpreter','latex','FontSize',18);
col=colorbar;
col.FontSize=16;
col.TickLabelInterpreter= 'latex';
% clim([0 1.1]*1e11)