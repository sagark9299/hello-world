function rho_l = density_hept(T)
C1 = 0.61259; C2 = 0.26211; C3 = 540.2; C4 = 0.28141; % constants for lambda0 calculation
rho_l = C1./C2.^(1+(1-(T./C3)).^C4); % in mol/dm^3
rho_l = rho_l.*1000;
end