clc; clear; close all;
nf = 500;
Lx0 = 0; Lx1 = 0.1; Ly0 = 0.1; Ly1 = 20; % 0.5 Lx1 = 0.1
dimtau = nf; dimxi = 50; %100
dtau = (Ly1-Ly0)/dimtau; dxi = (Lx1-Lx0)/dimxi;
tau = Ly0:dtau:Ly1;
xi = Lx0:dxi:Lx1;
xin = -1; xf = 10; ti = 0; tf = 2;
dimx = nf; dimt = 100;
dx = (xf-xin)/dimx; dt = (tf-ti)/dimt;
xv = xin:dx:xf; %0:1:20; 
D = 0.001;
V = 01;
%%
% f = 40;
% amp = 0.5;
% omega = 2*pi*f;
% [X,Y] = meshgrid(xi,tau);
% [X1,Y1] = meshgrid(xv,tv);
% L = length(xv);
% Id = zeros(1,L);
% phi = kaan_evap();
% for j = 1:1:length(tv)
%     for i = 1:L
%         t = tv(j);
%         x = xv(i);
%         G = heaviside(t-Y)./sqrt(pi*4*D.*(t-Y)).*(exp(-(X-x+V.*(t-Y)).^2./(4*D.*(t-Y))));
%         phi = exp(-(X-0.05).^2); %kaan_evap(tau) 
%         F = G.*phi;
%         I = trapz(tau,trapz(xi,F,2));
%         Id(i,j) = I;
%     end
% end
% Id_max = max(max(Id));
% plot(xv,Id(:,1)./Id_max,xv,Id(:,2)./Id_max)
% ylim([-1,1])
% grid on
%%
t = ti:dt:tf;
x = xin:dx:xf;
x1 = 0.1; x0 = 0;
xp = linspace(x0,x1,nf+1);
[X2,Y2] = meshgrid(xp,tau);
[X1,Y1] = meshgrid(x,tau);
% phi = square(X2)+1;
% G = 1./sqrt(pi*4*D.*Y1) .* (exp(-(X1-V.*Y1).^2./(4*D.*Y1))); 
G = -1./sqrt(pi*16*D.*Y1) .* (erf((X1-V.*Y1-x1)./(sqrt(4*D.*Y1)))-erf((X1-V.*Y1-x0)./(sqrt(4*D.*Y1))));%.* ...
   % exp(-0.1^2./(4*D.*Y1)); % earlier Y2 and X2
% G = -1./sqrt(pi*16*D.*tau) .* (erf((0.2-V.*tau-x1)./(sqrt(4*D.*tau)))-erf((0.2-V.*tau-x0)./(sqrt(4*D.*tau)))) .* ...
%     exp(-0.1^2./(4*D.*tau));
phi = kaan_evap(); 
F = G.*phi;
I = sum(F)*dtau;
% I = F.*dtau;
% figure(1); plot(tau,phi);
% figure(2)
% AxesH = axes('Xlim', [-1, 10], 'XTick', -1:1:20, 'NextPlot', 'add');
plot(x,I)
% set(gca,'XMinorGrid','on')
% set(gca,'YMinorGrid','on')
% ylim([0 3])
xlim([-1 xf])
grid on
% ax = gca;
% ax.FontSize = 15;
% ylabel ('$Y_{fuel}$','FontSize',15,'Interpreter','latex')
% xlabel('x [m]','FontSize',15,'Interpreter','latex')
% I(201)