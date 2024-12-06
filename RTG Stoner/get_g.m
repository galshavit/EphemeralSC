function [g,ind,phi_0] = get_g(F,phi,tuning_p)
minimas = 0*tuning_p; maximas = 0*tuning_p;
for i=1:length(tuning_p)
    f = F(i,:);
    [~,ind]=min(f);
    minimas(i) = phi(ind);
    maximas(i) = max(f(1:ind));
end
% jump = abs(diff(minimas));
% [~,ind] = max(jump);
% g = maximas(ind+1);
[g,indi]=max(maximas);
[~,I] = min(F(indi,:));
phi_0  = phi(I);
end