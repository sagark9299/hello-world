clear all; clc; close all;
L = 0.1;
n = 10;
T0 = 0;
T1s = 40;
T2s = 20;
dx = L/n; 
alpha = 0.0001; % thermal diffusivity
tf = 60;
dt = 0.1;
x = dx/2:dx:L-dx/2;
T = ones(n,1)*T0;
dTdt = zeros(n,1);
t = 0:dt:tf;
for i=1:length(t)
    for j = 2:n-1
        dTdt(j) = alpha*(-(T(j)-T(j-1))/dx^2+(T(j+1)-T(j))/dx^2);
    end
    dTdt(1) = alpha*(-(T(1)-T1s)/dx^2+(T(2)-T(1))/dx^2);
    dTdt(n) = alpha*(-(T(n)-T(n-1))/dx^2+(T2s-T(n))/dx^2);
    T = T+dTdt*dt;
    plot(x,T)
    pause(0.1)
end