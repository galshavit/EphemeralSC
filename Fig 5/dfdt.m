function dgdl = dfdt(f,l,tau,js,gamma)

dgdl = (-(js)^2 *exp(2*gamma*l)*4/(27*f^3) + f - f^3)/tau;
end