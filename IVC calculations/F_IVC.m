%%
clear
close all
clc
tic
%%
n = -1.5 *1e12;
Ds = linspace(23,25,181); % meV  V/nm --> D*3/200
g = 1.45*1e-11; % meV cm^2
numx=5800;
%
Delta = linspace(0,4,61);
temp = zeros(length(Ds),length(Delta));
mus=temp; Fkin=temp; clear('temp');
%
%%
for w=1:length(Ds)
    %%
    [E,x] = BandStructure6(numx,Ds(w));
    ek = squeeze(E(:,:,3)); clear('E');
    Omega = 1/(((x(1)-x(2))^2)/(4*pi^2))*1e4;
    %%
    for i=1:length(Delta)
        [mus(w,i),a,b] = get_mu_Delta(n,ek,Delta(i),Omega);
        Fkin(w,i) = (sum(a(a<=mus(w,i))) + sum(b(b<=mus(w,i))))/Omega ;
        disp(i);
        toc
    end
    disp(w);
    toc
end
save('FIVC')
%%
g = 1.46*1e-11; % meV cm^2
opt_D=0*Ds;
F = Fkin-Fkin(:,1)+repmat(Delta.^2,length(Ds),1)/g;
imagesc(Delta,Ds*3/200,F)
set(gca,'ydir','normal')
clim([-1 1]*2e9)
figure
hold on
for w=1:length(Ds)
    plot(Delta,Fkin(w,:)-Fkin(w,1)+Delta.^2/g,'.-')
    [~,ind]=min(Fkin(w,:)-Fkin(w,1)+Delta.^2/g);
    opt_D(w)=Delta(ind);
end