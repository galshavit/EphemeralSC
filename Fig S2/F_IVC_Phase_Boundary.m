%%
clear
close all
clc
tic
%%
n_list = linspace(-2*1e12,-1*1e12,26);
Ds = linspace(7,34,55); % meV  V/nm --> D*3/200
g = 1.46*1e-11; % meV cm^2
numx=1500;
%
Delta = linspace(0,4,51);
temp = zeros(length(n_list),length(Ds),length(Delta));
mus=temp; Fkin=temp; clear('temp');
%
%%
for w=1:length(Ds)
    %%
    [E,x] = BandStructure6(numx,Ds(w));
    ek = squeeze(E(:,:,3)); clear('E');
    Omega = 1/(((x(1)-x(2))^2)/(4*pi^2))*1e4;
    %%
    for z=1:length(n_list)
        n=n_list(z);
        parfor i=1:length(Delta)
            [mus(z,w,i),a,b] = get_mu_Delta(n,ek,Delta(i),Omega);
            Fkin(z,w,i) = (sum(a(a<=mus(z,w,i))) + sum(b(b<=mus(z,w,i))))/Omega ;
        end
    end
    disp(w);
    toc
end
save('ivc_boundary')
%%
opt_D = zeros(length(n_list),length(Ds));
F=0*Fkin;
for z=1:length(n_list)
    F(z,:,:) = squeeze(Fkin(z,:,:)-Fkin(z,:,1))+repmat(Delta.^2,length(Ds),1)/g;
end
for z=1:length(n_list)
    for w=1:length(Ds)
        [~,ind]=min(F(z,w,:));
        if ~isnan(mus(z,w,ind)); opt_D(z,w)=Delta(ind); end
    end
end

figure
imagesc(n_list,Ds*3/200,opt_D');
set(gca,'ydir','normal');
ylim([0.18 inf])
clim([1 inf])
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

