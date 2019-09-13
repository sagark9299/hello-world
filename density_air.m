function rho_0 = density_air(T)
mair = 28;
rho_0 = 360.77819.*(T.^1.00336); % in kg/m^3
rho_0 = rho_0*1000/mair; % in mol/m^3
end