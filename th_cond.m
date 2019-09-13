function kappa = th_cond(T)
c1k = 0.215; c2k = -0.000303;
kappa = c1k - c2k*T; %0.1171; %thermal conductivity of heptane at initial condition (323 K) % [W/m/K] Table 2-315 Perry
end