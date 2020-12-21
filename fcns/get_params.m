function p = get_params(gait)

p.predHorizon = 6;
p.simTimeStep = 1/200;
p.Tmpc = 4/100;             % MPC prediction step time
p.gait = gait;
p.Umax = 50;
p.decayRate = 1;
p.freq = 30;
p.Rground = eye(3);
p.Qf = diag([1e5 2e5 3e5 5e2 1e3 150 1e3 1e4 800 40 40 10]);

% ---- gait ----
if gait == 1                % 1 - bound
    p.Tst = 0.1;
    p.Tsw = 0.18;
    p.predHorizon = 7;
    p.simTimeStep = 1/100;  
    p.Tmpc = 2/100;
    p.decayRate = 1;
    p.R = diag(repmat([0.1 0.1 0.1]',[4,1]));
    p.Q = diag([5e4 2e4 1e6 4e3 5e2 5e2 1e4 5e4 1e3 1e2 5e2 1e2]);
    p.Qf = diag([2e5 5e4 5e6 8e3 5e2 5e2 1e4 5e4 5e3 1e2 1e2 1e2]);
elseif gait == 2            % 2 - pacing
    p.Tst = 0.12;
    p.Tsw = 0.12;
    p.R = diag(repmat([0.1 0.2 0.1]',[4,1]));
    p.Q = diag([5e3 5e3 9e4 5e2 5e2 5e2 7e3 7e3 7e3 5e1 5e1 5e1]);
elseif gait == 3            % 3 - gallop
    p.Tst = 0.08;
    p.Tsw = 0.2;
    p.R = diag(repmat([0.1 0.2 0.1]',[4,1]));
    p.Q = diag([3e3 3e3 4e6 5e2 1e3 150 1e4 1e4 800 1e2 5e1 5e1]);
elseif gait == 4            % 4 - trot run
    p.Tst = 0.12;
    p.Tsw = 0.2;
    p.Tmpc = 3/100;
    p.predHorizon = 6;
    p.decayRate = 1;
    p.R = diag(repmat([0.1 0.18 0.08]',[4,1]));
    p.Q = diag([1e5 1e5 1e5 1e3 1e3 1e3 2e3 1e4 800 100 40 10]);
    p.Qf = diag([1e5 1.5e5 2e4 1.5e3 1e3 100 2e3 2e3 800 100 60 10]);
elseif gait == 5            % 4 - crawl
    p.Tst = 0.3;
    p.Tsw = 0.1;
    p.R = diag(repmat([0.1 0.2 0.1]',[4,1]));
    p.Q = diag([5e5 5e5 9e5 5 5 5 3e3 3e3 3e3 3 3 3]);
else                        % 0 - trot
    p.predHorizon = 6;
    p.simTimeStep = 1/100;
    p.Tmpc = 8/100;
    p.Tst = 0.3;
    p.Tsw = 0.15;
    p.R = diag(repmat([0.1 0.2 0.1]',[4,1]));
    p.Q = diag([1e5 2e5 3e5 5e2 1e3 1e3 1e3 1e4 800 40 40 10]);
    p.Qf = p.Q;
end



%% Physical Parameters
p.mass = 5.5;
p.J = diag([0.026,0.112,0.075]);
p.g = 9.81;
p.mu = 1;       % friction coefficient
p.z0 = 0.2;     % nominal COM height
p.pf34 = [[0.15;0.094;0],[0.15;-0.094;0],[-0.15;0.094;0],[-0.15;-0.094;0]];

p.L = 0.301;    % body length
p.W = 0.088;    % body width
p.d = 0.05;     % ABAD offset
p.h = 0.05;     % body height
p.l1 = 0.14;    % link1 length
p.l2 = 0.14;    % link2 length

%% Swing phase 
p.Kp_sw = 300;  % Kp for swing phase

%% color
p.body_color    = [42 80 183]/255;
p.leg_color     = [7 179 128]/255;
p.ground_color  = [195 232 243]/255;



