function D = diffusion_coeff(T)
mw = 100.205; %molecular weight of heptane [g/mol]
mair = 28; %molecular weight of air [g/mol]
P = 1.01325; % bar, 1 atm
Mab = 2*((1/mw)+(1/mair))^-1;
sig_air = 3.62;
Vc = 428; % [cm^3/mol] critical volume of n-heptane
Vb = 0.285*(Vc)^1.048;
sig_hept = 1.18*Vb^(1/3);
sig_air_hept = (sig_air + sig_hept)/2;
ek_air = 97;
Tb = 371.6; % Normal boiling point of heptane
ek_hept = 1.15*Tb;
ek_air_hept = sqrt(ek_air*ek_hept);
Ts = T./ek_air_hept;
A = 1.06036; B = 0.15610; C = 0.19300; D = 0.47635; E = 1.03587; F = 1.52996; G = 1.76474; H = 3.89411;
omega = A./Ts.^B + C./exp(D.*Ts) + E./exp(F.*Ts) + G./exp(H.*Ts);
D = (3.03-(0.98./sqrt(Mab))).*1e-3.*T.^(3/2)./(P.*sqrt(Mab).*sig_air_hept.^2.*omega);
D = D*1e-4;
end