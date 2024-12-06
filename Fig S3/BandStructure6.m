function [E,x] = BandStructure6(numx,Delta)
%%
a = 2.46 *1e-10;
x=0.014*2*pi/a*linspace(-1,1,numx);
E = zeros(length(x),length(x),6);
%% meV
v0 = 3100 *sqrt(3)*a/2;
g1 = 380;
g2 = -15;
v3 = -290 *sqrt(3)*a/2;
v4 = -141 *sqrt(3)*a/2;
delta = -10.5;
Delta2 = -2.3;
%%
for i=1:length(x)
    for j=1:length(x)
        p = x(i)+1j*x(j);
        H=[Delta + Delta2 + delta, 0.5*g2, v0*conj(p), v4*conj(p),v3*p,0;
            0.5*g2, Delta2-Delta+delta, 0, v3*conj(p), v4*p, v0*p;
            v0*p , 0, Delta+Delta2, g1, v4*conj(p) ,0;
            v4*p, v3*p, g1, -2*Delta2, v0*conj(p), v4*conj(p);
            v3*conj(p), v4*conj(p), v4*p, v0*p, -2*Delta2, g1;
            0, v0*conj(p), 0, v4*p, g1, Delta2-Delta];
        E(i,j,:)=eig(H);
    end
end
end