%%
clear
close all
clc
%%
js=linspace(0.01,1.5,4001);
tau=1;
Opt    = odeset('Events', @myEvent,'relTol',1e-7);
gammas = [0.1 0.05 0.02 0.01 0.005];
td = 0*meshgrid(js,gammas);
td0=0*js;
%%
figure
hold on
for w=1:length(gammas)
    gamma=gammas(w);

    for i=1:length(js)
        [t,f] = ode23(@(t,f) dfdt(f,t,tau,js(i),gamma), [0,10000*tau], 1,Opt);
        td(w,i) = t(end);
    end
end
gamma=0;
js0=1+logspace(-5,-0.5,2001);
for i=1:length(js0)
    [t,f] = ode45(@(t,f) dfdt(f,t,tau,js0(i),gamma), [0,10000*tau], 1,Opt);
    td0(i) = t(end);
end
%%
figure
hold on
for w=1:length(gammas)
    tt=td(w,:);
    plot(js(js>sqrt(gammas(w))),tt(js>sqrt(gammas(w))),'.-','DisplayName', ['$\tau_{\rm decay}=\,$',sprintf('%d', 1/gammas(w)), '$\tau_{GL}$']  );
end
plot(js0(js0>1),td0(js0>1),'--k','LineWidth',2,'DisplayName','$\tau_{\rm decay}\to \infty$')
%
set(gca,'Yscale', 'log');
leg=legend;
leg.FontSize=14;
leg.Interpreter= 'latex';
box on
ax = gca;
ax.XAxis.FontSize = 22;
ax.YAxis.FontSize = 22;
set(gca,'TickLabelInterpreter', 'latex');
xlabel('$I/I_c^0$ ','interpreter','latex','FontSize',24);
ylabel('$t_{\rm delay}/\tau_{\rm GL}$','interpreter','latex','FontSize',24);
ylim([0 inf])
xlim([0 inf])
grid