function cp_0 = sp_heat_air(T)
C1 = 1.9327e-10; C2 = 7.9999e-7; C3 = 1.1407e-3; C4 = 4.4890e-1; C5 = 1.0575e3; % constants for Cp calculation
Mw = 28;
cp_0 = C1.*T^4 - C2.*T^3 + C3.*T^2 - C4.*T + C5; % J/mol*K Table 2-153 Perry Eq2 for heptane 234.12
cp_0 = cp_0.*Mw./1000;
end