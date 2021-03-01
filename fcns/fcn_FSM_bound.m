function [FSMout,Xd,Ud,Xt] = fcn_FSM_bound(t_,Xt,p)
%% parameters
acc = 1;
vd = [2;0];

Tst_ = p.Tst;
Tst = min(Tst_,0.2/norm(Xt(4:5)));
Tsw = p.Tsw;
T = Tst + Tsw;
Tair = 1/2 * (Tsw - Tst);

%%%%%%%% periodic traj for bounding %%%%%%%%%
p = fcn_bound_ref_traj(p);

%%%%%%%%% decompose state %%%%%%%%%%%%%%
[pc,dpc,vR,wb] = deal(Xt(1:3),Xt(4:6),Xt(7:15),Xt(16:18));
R = reshape(vR,[3,3]);
idx_pf = 19:30;
pf = Xt(idx_pf);
pf34 = reshape(pf,[3,4]);

%% initialization
persistent FSM Ta Tb pf_trans Ta_sw Tb_sw
if isempty(FSM)
    FSM = 1;
    Ta = 0;
    Tb = Ta + Tst;
    pf_trans = pf;
    Ta_sw(1:2) = Tst - T;
    Tb_sw(1:2) = Ta_sw(1:2) + Tsw;
    Ta_sw(3:4) = -Tair;
    Tb_sw(3:4) = Ta_sw(3:4) + Tsw;
end

%% FSM
% 1 - FrontStance -- Tst expires --> 2
% 2 - air_1 -- Tair expires --> 3
% 3 - BackStance -- Tst expires --> 4
% 4 - air_2 -- Tair expires --> 1

t = t_(1);
s = (t - Ta) / (Tb - Ta);
if (FSM == 1) && (s >= 1)
    Ta = Tb;
    Tb = Tb + Tair;
    FSM = FSM + 1;
    Ta_sw(1:2) = Ta_sw(1:2) + T;
    Tb_sw(1:2) = Ta_sw(1:2) + Tsw;
    pf_trans(1:6) = pf(1:6);
elseif (FSM == 2) && (s >= 1)
    Ta = Tb;
    Tb = Tb + Tst;
    FSM = FSM + 1;
elseif (FSM == 3) && (s >= 1)
    Ta = Tb;
    Tb = Tb + Tair;
    FSM = FSM + 1;
    Ta_sw(3:4) = Ta_sw(3:4) + T;
    Tb_sw(3:4) = Ta_sw(3:4) + Tsw;
    pf_trans(7:12) = pf(7:12);
elseif (FSM == 4) && (s >= 1)
    Ta = Tb;
    Tb = Tb + Tst;
    FSM = 1;
end

%% Xd/Ud
s = (t - Ta) / (Tb - Ta);
s(s<0) = 0;
s(s>1) = 1;

FSM_ = zeros(1,p.predHorizon);
s_ = zeros(1,p.predHorizon);

Ud = zeros(12,p.predHorizon);
Xd = repmat(Xt,[1,p.predHorizon]);

pd = [0;0];
dpd = [0;0];

for ii = 1:p.predHorizon
%%%%%%% Xd/Ud along the prediction horizon %%%%%%%%%%
[FSM_(ii),s_(ii)] = fcn_FSM_pred_hor(FSM,Ta,t_(ii),p);
FSM_(1) = FSM;
s_(1) = s;
for dir_xy = 1:2
    if t_(ii) < (vd(dir_xy) / acc)
        dpd(dir_xy) = acc * t_(ii);
        pd(dir_xy) = 1/2 * acc * t_(ii)^2;
    else
        dpd(dir_xy) = vd(dir_xy);
        pd(dir_xy) = vd(dir_xy) * t_(ii) - 1/2 * vd(dir_xy)^2/acc;
    end
end

if FSM_(ii) == 1                 % 1 - front stance
    Fz_d = polyval_bz(p.Fz_co,s_(ii));
    dz_d = polyval_bz(p.dz_co,s_(ii));
    z_d = polyval_bz(p.z_co,s_(ii));
    
    tau_d = polyval_bz(p.tau_co,s_(ii));
    dth_d = polyval_bz(p.dth_co,s_(ii));
    th_d = polyval_bz(p.th_co,s_(ii));
    R_d = expm(hatMap([0;th_d;0]));
    
    dx = pc(1) - pf34(1,1);
    Fx_d = (dx * Fz_d - tau_d) / (z_d);
    
    Ud(1:6,ii) = 1/2*[Fx_d;0;Fz_d;Fx_d;0;Fz_d];
    Xd(1:18,ii) = [pd;z_d;dpd;dz_d;R_d(:);0;dth_d;0];
elseif FSM_(ii) == 2             % 2 - air 1
    dz_d = p.dz_co(end) - p.g * (s_(ii) * Tair);
    z_d = p.z_co(end) + p.dz_co(end) * (s_(ii) * Tair) - ...
                        1/2 * p.g * (s_(ii) * Tair)^2;
    dth_d = p.dth_co(end);
    th_d = p.th_co(end) + dth_d * (s_(ii) * Tair);
    R_d = expm(hatMap([0;th_d;0]));
    Xd(1:18,ii) = [pd;z_d;dpd;dz_d;R_d(:);0;dth_d;0];
elseif FSM_(ii) == 3             % 3 - back stance
    Fz_d = polyval_bz(p.Fz_co,s_(ii));
    dz_d = polyval_bz(p.dz_co,s_(ii));
    z_d = polyval_bz(p.z_co,s_(ii));
    
    tau_d = polyval_bz(-p.tau_co,s_(ii));
    dth_d = polyval_bz(-p.dth_co,s_(ii));
    th_d = polyval_bz(-p.th_co,s_(ii));
    R_d = expm(hatMap([0;th_d;0]));
    
    dx = pc(1) - pf34(1,3);
    Fx_d = (dx * Fz_d - tau_d) / (z_d);
    
    Ud(7:12,ii) = 1/2*[Fx_d;0;Fz_d;Fx_d;0;Fz_d];
    Xd(1:18,ii) = [pd;z_d;dpd;dz_d;R_d(:);0;dth_d;0];
elseif FSM_(ii) == 4             % 4 - air 2
    dz_d = p.dz_co(end) - p.g * (s_(ii) * Tair);
    z_d = p.z_co(end) + p.dz_co(end) * (s_(ii) * Tair) - ...
                        1/2 * p.g * (s_(ii) * Tair)^2;
    dth_d = -p.dth_co(end);
    th_d = -p.th_co(end) + dth_d * (s_(ii) * Tair);
    R_d = expm(hatMap([0;th_d;0]));
    Xd(1:18,ii) = [pd;z_d;dpd;dz_d;R_d(:);0;dth_d;0];
end

end


%% swing leg
[L,W,d] = deal(p.L,p.W,p.d);
p_hip_b = [[L/2;W/2+d;0],[L/2;-W/2-d;0],[-L/2;W/2+d;0],[-L/2;-W/2-d;0]];
p_hip_R = R * p_hip_b;
ws = R * wb;
v_hip_R = repmat(dpc,[1,4]) + hatMap(ws) * p_hip_R;

% capture point
% p_cap = Tsw/2 * vd + sqrt(z0/g) * (vc - vd)
p_cap = zeros(2,4);
dpd = Xd(4:5,1);
for i_leg = 1:4
    temp = 0.8 * Tst * dpd + sqrt(p.z0/p.g) * (v_hip_R(1:2,i_leg) - dpd);
    temp(temp < -0.15) = -0.15;
    temp(temp > 0.15) = 0.15;
    p_cap(:,i_leg) = pc(1:2) + p_hip_R(1:2,i_leg) + temp;
end

% desired foot placement
pfd = Xt(idx_pf);
s_sw = (t - Ta_sw) ./ (Tb_sw - Ta_sw);
s_sw(s_sw<0) = 0;s_sw(s_sw>1) = 1;

if FSM == 1
    for i_leg = 3:4
        idx = 3*(i_leg-1) + (1:3);
        co_x = linspace(pf_trans(idx(1)),p_cap(1,i_leg),6);
        co_y = linspace(pf_trans(idx(2)),p_cap(2,i_leg),6);
        co_z = [0 0 0.15 0.15 0 -0.002];
        pfd(idx) = [polyval_bz(co_x,s_sw(i_leg));
                    polyval_bz(co_y,s_sw(i_leg));
                    polyval_bz(co_z,s_sw(i_leg))];
    end
elseif FSM == 2
    for i_leg = 1:4
        idx = 3*(i_leg-1) + (1:3);
        co_x = linspace(pf_trans(idx(1)),p_cap(1,i_leg),6);
        co_y = linspace(pf_trans(idx(2)),p_cap(2,i_leg),6);
        co_z = [0 0 0.15 0.15 0 -0.002];
        pfd(idx) = [polyval_bz(co_x,s_sw(i_leg));
                    polyval_bz(co_y,s_sw(i_leg));
                    polyval_bz(co_z,s_sw(i_leg))];
    end
elseif FSM == 3
    for i_leg = 1:2
        idx = 3*(i_leg-1) + (1:3);
        co_x = linspace(pf_trans(idx(1)),p_cap(1,i_leg),6);
        co_y = linspace(pf_trans(idx(2)),p_cap(2,i_leg),6);
        co_z = [0 0 0.15 0.15 0 -0.002];
        pfd(idx) = [polyval_bz(co_x,s_sw(i_leg));
                    polyval_bz(co_y,s_sw(i_leg));
                    polyval_bz(co_z,s_sw(i_leg))];
    end
elseif FSM == 4
    for i_leg = 1:4
        idx = 3*(i_leg-1) + (1:3);
        co_x = linspace(pf_trans(idx(1)),p_cap(1,i_leg),6);
        co_y = linspace(pf_trans(idx(2)),p_cap(2,i_leg),6);
        co_z = [0 0 0.15 0.15 0 -0.002];
        pfd(idx) = [polyval_bz(co_x,s_sw(i_leg));
                    polyval_bz(co_y,s_sw(i_leg));
                    polyval_bz(co_z,s_sw(i_leg))];
    end
end

Xd(idx_pf,:) = repmat(pfd,[1,p.predHorizon]);
Xd = repmat(Xd(:,1),[1,p.predHorizon]);

%% output
FSMout = FSM;

end

function [FSMout,sout] = fcn_FSM_pred_hor(FSM,Ta,t,p)
Tst = p.Tst;
Tsw = p.Tsw;
T = Tst + Tsw;
Tair = 1/2 * (Tsw - Tst);

tp = mod(t - Ta,T);

if (FSM == 1) || (FSM == 3)
    Tnode = [Tst Tair Tst Tair];
else
    Tnode = [Tair Tst Tair Tst];
end

for ii = 1:4
    if tp <= sum(Tnode(1:ii))
        FSMout = mod(FSM + ii - 1,4);
        sout = (tp - sum(Tnode(1:ii-1))) / Tnode(ii);
        break;
    end
end

end




