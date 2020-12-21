function [p,Xt,Ut] = fcn_bound_ref_traj(p)
% This function finds the initial condition for periodic bounding
% The calculation is based on the paper (citation):
% Park, Hae-Won, Patrick M. Wensing, and Sangbae Kim. 
% "High-speed bounding with the MIT Cheetah 2: Control design and experiments."
% The International Journal of Robotics Research 36, no. 2 (2017): 167-192.

[mass,J,g,Tst,Tsw] = deal(p.mass,p.J,p.g,p.Tst,p.Tsw);
T = Tst + Tsw;
Tair = 1/2 * (Tsw - Tst);

b_co = [0 0.8 1 1 0.8 0];
b_ = mean(b_co);

%% Fz
% 2 * alpha * b_ * Tst = mass * g * T
alpha_z = (mass * g * T) / (2 * b_ * Tst);
Fz_co = alpha_z * b_co;

dz_co = bz_int(Fz_co/mass-g,0,Tst);
z_co = bz_int(dz_co,0,Tst);

% first principle: integration
dz0 = -1/(Tst+Tair)*(z_co(end) + Tair*(dz_co(end)+g*Tst)-1/2*g*((Tst+Tair)^2-Tst^2));

dz_co = bz_int(Fz_co/mass-g,dz0,Tst);
z_co = bz_int(dz_co,p.z0,Tst);

%% theta
alpha_th = 140 * J(2,2);
tau_co = -alpha_th * b_co;

dth_co = bz_int(tau_co/J(2,2),0,Tst);
% by symmetry
dth0 = -1/2 * dth_co(end);

th0 = dth0*Tair/2;
dth_co = bz_int(tau_co/J(2,2),dth0,Tst);
th_co = bz_int(dth_co,th0,Tst);

%% output B-spline coefficient
p.Fz_co = Fz_co;
p.dz_co = dz_co;
p.z_co = z_co;

p.tau_co = tau_co;
p.dth_co = dth_co;
p.th_co = th_co;

%% intial condition
R0 = expm(hatMap([0;th0;0]));
Xt = [0;0;p.z0;0;0;dz0;R0(:);0;dth0;0];
Xt(19:30) = p.pf34(:);
Ut = repmat([0;0;1/4*p.mass*p.g],[4,1]);

end