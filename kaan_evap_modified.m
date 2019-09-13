% function mass = kaan_evap()
clear ;clc; close all;
t = linspace(0,15,51);
x = linspace(0,0.1,51);
[X,Y] = meshgrid(x,t);
% t = T;
% t = 0:0.05:5;
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
D0 = 100e-6; % initial liquid droplet size % [m]
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
Tl =  Tin - C2 - ((Tin-Tl0)-C2)*exp(-C1.*Y/th);
% Tl = Tl0 - (C1*C2.*t)./th;
Twb = Tin - (BT0*lambda0)/cp_0; % wet bulb temperature (i.e. Tl when t goes to infinity)
tau_wb = (th/C1)*log((Tin-Tl0)-C2/(Tin-Twb)-C2);
%%
% uc= 0.5;
% U=10;
% up0=10;
% T_vis = 273.15;
% muc_ref = 1.716e-5;
% S = 110.4;                                 
% muc = muc_ref*((Tref/T_vis)^(3/2))*((T_vis+S)/(Tref+S));
rho_air = 0.0353e3; % density of air at Tref in mol/m^3 0.0342e3 0.0353e3(313)
% nuc = muc/(rho_air*mair/1000);
% Rep_max = (D0*abs(U+uc+up0))/(nuc);
p_ratio_tr = exp((h_vap(Tl)./R).*((1./Tb)-(1./Tl)));
Xs_tr = p_ratio_tr;
mmix_tr =  mw.*Xs_tr + mair.*(1-Xs_tr);
Ys_tr = Xs_tr.*mw./mmix_tr;
BM_tr = (Ys_tr-Yin)./(1-Ys_tr);
dt = t(2)-t(1);
r2 = zeros(size(t,2));
%%
% r0 = ones(1,size(t,2))*D0/2;
r0 = D0;
r2 = ones(size(t,2))*r0^2;
K_tr = -8.*rho_air.*1.4e-8.*log(1+BM_tr)./density_hept(Tl); %diffusion_coeff(Tref_tr) 1.25e-8
% K_conv = K_tr * (1+0.3*Rep_max^0.28);
for j=1:length(t)
    for i=1:length(t)-1
    r2(i+1,j) = r2(i) + K_tr(i)*dt;
    end
end
idx = r2 <= 0;
r2(idx) = 0;
% rn=r2./r0^2;
% I = find(abs(t-tau_wb)<1e-3);
r = sqrt(r2);
% figure(2)
% plot(t,rn)
% ylim([0,1])
%%
rhod=980;
rhoc=1.2;
gamma=rhod/rhoc;
muc=1.8559e-05;
nuc=muc/rhoc;
f=250; 
omega=2*pi*f;
uc= 0.5;
U=2;
up0=2;
syms n

% x = linspace(0,0.1,51);
% x = X;
% [X,Y] = meshgrid(x,t);
C=18*nuc./(gamma.*r.^2);
Uc = U;
up = up0;
% for j = 1:length(t)
%     for i = 1:length(x)
%         C = 18*nuc/(gamma*r(j)^2);
%        
%         z(i) = exp(-(Uc-up+(C*x(i)))/Uc)*((up/Uc) -1);
        z = exp(-(Uc-up+(C.*X))/Uc).*((up/Uc) -1);

%         W(i) = symsum(((((-n)^(n-1))/factorial(n))*(z(i)^n)),n,1,50);
        W = symsum(((((-n)^(n-1))/factorial(n)).*(z.^n)),n,1,50);
% 
%         ts(i) = (x(i)/Uc)+ ((Uc-up)/(C*Uc))+ (W(i)/C);
        ts = (X./Uc)+ ((Uc-up)./(C.*Uc))+ (W./C);
% 
%         A(i) = -(Uc-up)*exp(-C*ts(i));
        A = -(Uc-up).*exp(-C.*ts);
%         
        B = Uc ;
% 
%         D(i) = uc*(1-exp(-C*ts(i)));
        D = uc.*(1-exp(-C.*ts));
% 
%         k = 1;   %1
%       %%ask !! omega in paranthesis
%         E(i) = (((omega/(gamma*C))*cos((omega*ts(i))-(2*pi)))- sin((omega*ts(i))-(2*pi)));
        E = (((omega./(gamma.*C)).*cos((omega.*ts)-(2*pi)))- sin((omega.*ts)-(2*pi)));
% 
%         ND1(i) = 1./(A(i) + B + (D(i)*E(i)));
        ND1 = 1./(A + B + (D.*E));
% 
%         ND2(i) = 1/(A(i) + B);
        ND2 = 1./(A + B);
% 
%         ND(i) = ND1(i)/ND2(i) ;
        ND = ND1./ND2 ;
%         IND(i,j) = ND(i);
%     end
% end
mass = double(ND);
% figure(1)
% plot(x,ND,'MarkerSize',4);
% ax = gca;
% ax.FontSize = 15;
% ylabel ('$\overline{ND}$','FontSize',15,'Interpreter','latex')
% xlabel('x [m]','FontSize',15,'Interpreter','latex')
% xlim([0 0.04])
%%
% mdot = pi.*980.*5.5e-8.*ND.*r/4;
mdot = -pi.*rho_l0_mass.*K_tr(end).*ND.*r/4;
mdot_norm = mdot./max(mdot);
mliq = double(mdot_norm);
% mvap = (1 - mliq).*double(ND);
% mass = mliq;
% plot(x,mdot*1000)
% grid on
% ax = gca;
% ax.FontSize = 15;
% ylabel ('$\dot{m}[kg/s]$','FontSize',15,'Interpreter','latex')
% xlabel('x [m]','FontSize',15,'Interpreter','latex')
% xlim([0 0.08])
% end