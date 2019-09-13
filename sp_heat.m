function cp_l = sp_heat(T)
C1 = 61.26; C2 = 314410; C3 = 1824.6; C4 = -2547.9; % constants for Cp calculation
Tr = T./540.2;
t = 1-Tr;
cp_l = C1^2/t + C2 - 2*C1*C3*t - C1*C4*t^2 - C3^2*t^3/3 - C3*C4*t^4/2 - C4^2*t^5/5; % J/mol*K Table 2-153 Perry Eq2 for heptane 234.12
cp_l = cp_l/1000;
end