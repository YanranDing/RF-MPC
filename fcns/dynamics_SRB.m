function dXdt = dynamics_SRB(t,Xt,Ut,Xd,U_ext,p)

%% parameters
mass = p.mass;
J = p.J;       % inertia tensor in body frame {B}
g = 9.81;

%% decompose
% X = [pc dpc vR wb pf]'
pc = reshape(Xt(1:3),[3,1]);
dpc = reshape(Xt(4:6),[3,1]);
R = reshape(Xt(7:15),[3,3]);
wb = reshape(Xt(16:18),[3,1]);
pf34 = reshape(Xt(19:30),[3,4]);
pfd34 = reshape(Xd(19:30),[3,4]);

% r
r34 = pf34 - repmat(pc,[1,4]);

% GRF
f34 = reshape(Ut,[3,4]);

%% dynamics
ddpc = 1/mass * (sum(f34,2) + U_ext) + [0;0;-g];
dR = R * hatMap(wb);

tau_s = zeros(3,1);       % body torque expressed in {S}
for ii = 1:4
    tau_s = tau_s + hatMap(r34(:,ii)) * f34(:,ii);
end
tau_ext = hatMap(R * p.p_ext) * U_ext;
tau_tot = sum(tau_s,2) + tau_ext;
dwb = J \ (R' * tau_tot - hatMap(wb) * J * wb);

dpf = p.Kp_sw * (pfd34(:) - pf34(:));

dXdt = [dpc;ddpc;dR(:);dwb;dpf];


end







