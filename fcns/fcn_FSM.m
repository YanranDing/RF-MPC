function [FSMout,Xd,Ud,Xt] = fcn_FSM(t_,Xt,p)

%% parameters
[L,W,d] = deal(p.L,p.W,p.d);
gait = p.gait;
Tst_ = p.Tst;
Tst = min(Tst_,0.2/norm(Xt(4:5)));
Tsw = p.Tsw;
T = Tst + Tsw;
Tair = 1/2 * (Tsw - Tst);

[pc,dpc,vR,wb] = deal(Xt(1:3),Xt(4:6),Xt(7:15),Xt(16:18));
R = reshape(vR,[3,3]);
idx_pf = 19:30;
pf34 = reshape(Xt(idx_pf),[3,4]);

%% initialization
persistent FSM Ta Tb pf_R_trans
if isempty(FSM)
    FSM = zeros(4,1);
    Ta = zeros(4,1);
    Tb = ones(4,1);
    pf_R_trans = Xt(idx_pf);
end

t = t_(1);      % current time
s = zeros(4,1);

%% FSM
% 1 - stance
% 2 - swing

for i_leg = 1:4
    s(i_leg) = (t - Ta(i_leg)) ./ (Tb(i_leg) - Ta(i_leg));
    s(s<0) = 0;
    s(s>1) = 1;
    % --- FSM ---
    if FSM(i_leg) == 0          % init to stance
        if gait == -3           % backflip
            Ta(i_leg) = t;
            Tb([1,2]) = Ta([1,2]) + p.Tds*[1;1];
            Tb([3,4]) = Ta([3,4]) + p.Tbs*[1;1];
        elseif gait == -1           % pose control
            Ta(i_leg) = 0;
            Tb(i_leg) = -1;
        elseif gait == 1        % bound
            Ta([1,2]) = [t;t];
            Ta([3,4]) = [1;1] * (t + 1/2*(Tst + Tsw));
            Tb([1,2]) = Ta([1,2]) + Tst;
            Tb([3,4]) = Ta([3,4]) + Tst + Tair;
        elseif gait == 2        % pacing
            Ta([1,3]) = [t;t];
            Ta([2,4]) = [1;1] * (t + 1/2*(Tst + Tsw));
            Tb(i_leg) = Ta(i_leg) + Tst;
        elseif gait == 3        % gallop
            Ta(1) = t;
            Ta(2) = t + 0.05;
            Ta(3) = t + 0.05 + Tst;
            Ta(4) = t + 0.1 + Tst;
            Tb(i_leg) = Ta(i_leg) + Tst;
        elseif gait == 5        % crawl
            Ta(1) = t;
            Ta(2) = t + Tsw;
            Ta(3) = t + Tsw*2;
            Ta(4) = t + Tsw*3;
            Tb(i_leg) = Ta(i_leg) + Tst;
        else                    % trot walk
            Ta([1,4]) = [t;t];
            Ta([2,3]) = [1;1] * (t + 1/2*(Tst + Tsw));
            Tb(i_leg) = Ta(i_leg) + Tst;
        end
        FSM(i_leg) = FSM(i_leg) + 1;
        pf_R_trans = Xt(idx_pf);
    elseif FSM(i_leg) == 1 && (s(i_leg) >= 1 - 1e-7)    % stance to swing
        FSM(i_leg) = FSM(i_leg) + 1;
        Ta(i_leg) = t;
        Tb(i_leg) = Ta(i_leg) + Tsw;
        pf_R_trans = Xt(idx_pf);
    elseif FSM(i_leg) == 2 && (s(i_leg) >= 1 - 1e-7)    % swing to stance
        if gait == -3
            ph = R * [[L/2;W/2+d;0],[L/2;-W/2-d;0]] + repmat(pc,[1,2]);
            pf = ph + repmat([0;0;-0.18],[1,2]);
            if pf(3,1) < 1e-4
                FSM(i_leg) = 1;
                Ta(i_leg) = t;
                Tb(i_leg) = Ta(i_leg) + Tst;
                pf_R_trans = Xt(idx_pf);
            end
        else
            FSM(i_leg) = 1;
            Ta(i_leg) = t;
            Tb(i_leg) = Ta(i_leg) + Tst;
            pf_R_trans = Xt(idx_pf);
        end
    end
end

s = (t - Ta) ./ (Tb - Ta);
s(s<0) = 0;
s(s>1) = 1;

%% FSM in prediction horizon
FSM_ = repmat(FSM,[1,p.predHorizon]);
for i_leg = 1:4
    for ii = 2:p.predHorizon
        if t_(ii) <= Ta(i_leg)
            FSM_(i_leg,ii) = 1;
        elseif (Ta(i_leg) < t_(ii)) && (t_(ii) < Tb(i_leg))
            FSM_(i_leg,ii) = FSM(i_leg);
        elseif Ta(i_leg) + Tst + Tsw < t_(ii)
            FSM_(i_leg,ii) = FSM(i_leg);
        else
            if FSM(i_leg) == 1
                FSM_(i_leg,ii) = 2;
            else
                FSM_(i_leg,ii) = 1;
            end
        end
    end
end

if gait == -1       % pose
    FSM_ = ones(size(FSM_));
elseif gait == -3   % backflip
    FSM(3:4) = 1;
    FSM_(:,3:4) = 1;
end

% [4,predHorizon]: bool matrix
bool_inStance = (FSM_ == 1);

%% Gen ref traj
% --- Xd/Ud ---
[Xd,Ud] = fcn_gen_XdUd(t_,Xt,bool_inStance,p);

p = fcn_bound_ref_traj(p);

if gait == 1        % bound
    for ii = 1:p.predHorizon
        fsm = FSM_(:,ii);
        if fsm(1) == 1      % front stance
            s_ph = (t_(ii) - Ta(1)) / (Tb(1) - Ta(1));
            
            th_d = polyval_bz(-p.th_co,s_ph);
            dth_d = polyval_bz(-p.dth_co,s_ph);
            z_d = polyval_bz(p.z_co,s_ph);
            vR_d = reshape(expm(hatMap([0;th_d;0])),[9,1]);
            Xd(3,ii) = z_d;
            Xd(7:15,ii) = vR_d;
            Xd(17,ii) = dth_d;
            
            Fz_d = polyval_bz(p.Fz_co,s_ph);
            tau_d = polyval_bz(p.tau_co,s_ph);
            r = pf34(:,1) - pc;
            Ud([3,6],ii) = 1/2*Fz_d;
            Ud([1,4],ii) = 1/2*(r(1)*Fz_d - tau_d) / r(3);
        elseif fsm(3) == 1  % back stance
            s_ph = (t_(ii) - Ta(3)) / (Tb(3) - Ta(3));
            
            th_d = polyval_bz(p.th_co,s_ph);
            dth_d = polyval_bz(p.dth_co,s_ph);
            z_d = polyval_bz(p.z_co,s_ph);
            vR_d = reshape(expm(hatMap([0;th_d;0])),[9,1]);
            Xd(3,ii) = z_d;
            Xd(7:15,ii) = vR_d;
            Xd(17,ii) = dth_d;
            
            Fz_d = polyval_bz(p.Fz_co,s_ph);
            tau_d = polyval_bz(-p.tau_co,s_ph);
            r = pf34(:,3) - pc;
            Ud([9,12],ii) = 1/2*Fz_d;
            Ud([7,10],ii) = 1/2*(r(1)*Fz_d - tau_d) / r(3);
        end
    end
end

%% swing leg
p_hip_b = [[L/2;W/2+d;0],[L/2;-W/2-d;0],[-L/2;W/2+d;0],[-L/2;-W/2-d;0]];
p_hip_R = R * p_hip_b;
ws = R * wb;
v_hip_R = repmat(dpc,[1,4]) + hatMap(ws) * p_hip_R;

% capture point
% p_cap = Tsw/2 * vd + sqrt(z0/g) * (vc - vd)
p_cap = zeros(2,4);
vd = Xd(4:5,1);
for i_leg = 1:4
    temp = 0.8 * Tst * vd + sqrt(p.z0/p.g) * (v_hip_R(1:2,i_leg) - vd);
    temp(temp < -0.15) = -0.15;
    temp(temp > 0.15) = 0.15;
    p_cap(:,i_leg) = pc(1:2) + p_hip_R(1:2,i_leg) + temp;
end

% desired foot placement
if p.gait == -3     % backflip
    pfd = Xt(idx_pf);
    for i_leg = 1:2
        if FSM(i_leg) == 2
            ph = R * [[L/2;W/2+d;0],[L/2;-W/2-d;0]] + repmat(pc,[1,2]);
            pf = ph + repmat([0;0;-0.18],[1,2]);
            pf(3,:) = max([0,0],pf(3,:));
            pfd(1:6) = pf(:);
        end
    end
    Xd(idx_pf,:) = repmat(pfd,[1,p.predHorizon]);
elseif p.gait == -2     % GOT
    Rg = fcn_GOT_Rg_BB(t,p);
    Xt(idx_pf) = reshape(Rg * p.pf34,[12,1]);
    Xd(idx_pf,:) = repmat(Xt(idx_pf),[1,p.predHorizon]);
else
    pfd = Xt(idx_pf);
    for i_leg = 1:4
        idx = 3*(i_leg-1) + (1:3);
        if FSM(i_leg) == 2
            co_x = linspace(pf_R_trans(idx(1)),p_cap(1,i_leg),6);
            co_y = linspace(pf_R_trans(idx(2)),p_cap(2,i_leg),6);
            co_z = [0 0 0.1 0.1 0 -0.002];
            pfd(idx) = [polyval_bz(co_x,s(i_leg));
                        polyval_bz(co_y,s(i_leg));
                        polyval_bz(co_z,s(i_leg))];
        end
    end
    Xd(idx_pf,:) = repmat(pfd,[1,p.predHorizon]);
end

%% output
FSMout = FSM;
end

















