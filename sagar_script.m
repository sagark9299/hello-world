% function S_evap = sagar_script()
clear; clc; close;
t = linspace(0,10,5001);
Tb = 371.6; % Normal boiling point of heptane
mw = 100.205; %molecular weight of heptane [g/mol]
mair = 28; %molecular weight of air [g/mol]
Tin = 473.15; % gas temperature at infinity
Tl0 = 313.15; % initial liquid temperature
Tc = 540.2; % critical temperature of heptane
Tr = Tl0/Tc; % reduced temperature
Tref = Tl0+(1/3)*(Tin-Tl0); % reference temperature for gas properties calculation
kap0 = th_cond_air(Tref); %thermal conductivity of air at Tref % [W/m/K] (obtained from interpolating between 300 & 500 K at 0.1 MPa, Table 2-187) 31.343e-3 30.439e-3 (313)
kap_l0 = th_cond(Tl0);
cp_l0 =  sp_heat(Tl0); %234.12; % J/mol*K Table 2-153 Perry Eq2 for heptane 234.12 sp_heat_air(Tref)
cp_l0_mass = cp_l0/mw; % kJ/(kg*K)
cp_0 = sp_heat_air(Tref)+1; % J/mol*K at Tref (obtained from interpolating between 300 & 500 K at 0.1 MPa, Table 2-187) 0.029351e3 0.029352e3
cp_0_mass = cp_0/mair; % kJ/(kg*K)
rho_l0 = density_hept(Tl0); %liquid density at 323 K [mol/m^3]
rho_l0_mass = rho_l0*mw/1000; % [kg/m^3]
D_AB = diffusion_coeff(Tref); % diffusion coefficient of heptane in air [m^2/s] at Tref from Wilke-Lee model
% D_AB = 8e-1;
alpha_l0 = kap_l0/(rho_l0*cp_l0); % thermal diffusivity at 323 K %[m^2/s]
Le0 = 2.8976; % 2.8976
lambda0 = h_vap(Tl0); % heat of vaporization at 323 K J/mol
lambda0_mass = lambda0/mw; % [kJ/kg]
R = 8.314; % J/mol K
D0 = 20e-6; % initial liquid droplet size % [m]
th = D0^2/(4*alpha_l0);
p_ratio = exp((lambda0/R)*((1/Tb)-(1/Tl0))); % Clausius - Clapeyron: p/p_atm
Ys0_mol = p_ratio; % mole fraction - isobaric evaporation
mmix = mw*Ys0_mol + mair*(1-Ys0_mol); % mixture molecular weight
Ys_0_mass = Ys0_mol*mw/mmix;
Yin = 0; %mass fraction of fuel at infinity
Y_ref_mass = Ys_0_mass + (1/3)*(Yin - Ys_0_mass);
BM0 = Ys_0_mass-Yin/(1-Ys_0_mass);
BT0 = (1+BM0)^(1/Le0)-1;
C2 = lambda0*BT0/cp_0;
C1 = log(1+BM0)*kap0/(Le0*BT0*kap_l0);
Tl =  Tin - C2 - ((Tin-Tl0)-C2)*exp(-C1.*t/th);
% Tl = Tl0 - (C1*C2.*t)./th;
Twb = Tin - (BT0*lambda0)/cp_0; % wet bulb temperature (i.e. Tl when t goes to infinity)
tau_wb = (th/C1)*log((Tin-Tl0)-C2/(Tin-Twb)-C2);
% I = find(abs(Tl-Twb)<=1e-6);
% I1 = I;
% t = t(1:I1);
% tau_wb = th*(Tl0-Twb)/(C1*C2);
% figure(1)
% plot(t,Tl)
% grid on
% ylabel ('T_{Liq}','FontSize',15)
% xlabel('t [s]','FontSize',15)
% ylim([300 325])
%%
p_ratio_tr = exp((h_vap(Tl)/R).*((1./Tb)-(1./Tl)));
Xs_tr = p_ratio_tr;
mmix_tr =  mw.*Xs_tr + mair.*(1-Xs_tr);
Ys_tr = Xs_tr.*mw./mmix_tr;
BM_tr = (Ys_tr-Yin)./(1-Ys_tr);
dt = t(2)-t(1);
r2 = zeros(1,size(t,2));
r0 = D0;
r2(1) = r0^2;
Tref_tr = Tl+(1/3).*(Tin-Tl); 
rho_air = 0.0353e3; % density of air at Tref in mol/m^3 0.0342e3 0.0353e3(313)
K_tr = -8.*rho_air.*1.4e-8.*log(1+BM_tr)./density_hept(Tl); %diffusion_coeff(Tref_tr) 1.25e-8
% I = find (diff(K_tr)==0);
for i=1:length(t)-1
    r2(i+1) = r2(i) + K_tr(i)*dt;
%     if r2(i+1) < 0
%         break
%     end
end
idx = r2 <= 0;
r2(idx) = 0;
rn=r2./r0^2;
r = sqrt(r2);
% plot(t,rn)
% ylim([0,1])
% % xlim([0 10])
% grid on;
plot(t,rn)
grid on
ax = gca;
ax.FontSize = 15;
ylabel ('K (m^2/s)','FontSize',15)
xlabel('t [s]','FontSize',15)
%% after reaching Twb
% p_ratio_twb = exp((h_vap(Twb)/R)*((1/Tb)-(1/Twb)));
% Xs_twb = p_ratio_twb;
% mmix_twb =  mw*Xs_twb + mair*(1-Xs_twb);
% Ys_twb = Xs_twb*mw/mmix_twb;
% BM_twb = (Ys_twb-Yin)/(1-Ys_twb);
% Tref_wb = Twb+(1/3).*(Tin-Twb); % reference temperature for gas properties calculation
% K_twb =  2*rho_air*0.25e-7*log(1+BM_twb)/density_hept(Twb);
% t1 = t(I:end);
% r2a = zeros(1,size(t1,2));
% r2a(1) = r2(I(1));
% for k=1:length(t1)-1
%     r2a(k+1) = r2a(k) - K_twb*dt;
% %     if r2a(j+1) <= 0
% %         break
% %     end
% end
% % r2a = r2(I)^2 - K_twb*t1;
% rn1 = r2a./r0^2;
% % figure(3)
% plot(t1,rn1,t(1:I),rn(1:I),'b')
% ylim([0 1.2])
% xlim([0 10])
% grid on
% ylabel ('(d/d_{0})^2')
% xlabel('t [s]')
% tlife = tau_wb + r2a(1)/K_twb;
% S_evap = K_tr;
% end