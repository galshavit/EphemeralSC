%%
clear
close all
clc
load('BBG_DOS.mat')
tic
%%
ntot = -6.2 *1e11;
d=2;
delta = linspace(0,abs(ntot)/(2*d),3000);
%%
Fkin = zeros(length(Ds),length(delta));
Delta = repmat(delta,length(Ds),1);
%%
for i=1:length(Ds)
    n_eps = n_of_eps(eps,dos(i,:));
    for j=1:length(delta)
        Fkin(i,j)=Generic_Stoner_Fkin(eps,dos(i,:),n_eps + cumsum(n_eps*0+0.01),ntot,d,delta(j));
    end
    disp(i);
    toc
end
%%
U = 3.6e-11;
F = Fkin - U*Delta.^2;
figure
imagesc(delta,Ds,F);
clim([-2 2]*1e9)
colormap(jet)