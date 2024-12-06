function [n_eps] = n_of_eps(eps,dos)
% Relative to CNP, valence band version (so all negatives)
de = eps(2)-eps(1);
n_eps = de*cumsum(dos) - de*sum(dos);
end