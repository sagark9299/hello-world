function lambda = h_vap(T)
c1l = 5.0014e+7; c2l = 0.38795; % constants for lambda0 calculation
Tr = T./540.2;
lambda = (c1l.*(1-Tr).^c2l)./1000; % heat of vaporization at 323 K J/mol
end