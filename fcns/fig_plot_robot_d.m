function fig_plot_robot_d(Xd,Ud,p)

%% parameters
L = p.L;
W = p.W;
h = p.h;

body_color = p.body_color;
leg_color = p.leg_color;
ground_color = p.ground_color;

%% unpack
% X = [pc dpc vR wb pf]'
pcom =    reshape(Xd(1:3),[3,1]);
dpc =   reshape(Xd(4:6),[3,1]);
R =     reshape(Xd(7:15),[3,3]);
wb =    reshape(Xd(16:18),[3,1]);
pf34 =  reshape(Xd(19:30),[3,4]);

% GRF
f34 = reshape(Ud,[3,4]);


%% forward kinematics
% hips
Twd2com = [R,    pcom;
           0 0 0 1];
Tcom2h1 = [eye(3) [L/2 W/2 0]';
            0 0 0 1];
Tcom2h2 = [eye(3) [L/2 -W/2 0]';
            0 0 0 1];
Tcom2h3 = [eye(3) [-L/2 W/2 0]';
            0 0 0 1];
Tcom2h4 = [eye(3) [-L/2 -W/2 0]';
            0 0 0 1];
Twd2h1 = Twd2com * Tcom2h1;
Twd2h2 = Twd2com * Tcom2h2;
Twd2h3 = Twd2com * Tcom2h3;
Twd2h4 = Twd2com * Tcom2h4;

p_h1_wd = Twd2h1(1:3,4);
p_h2_wd = Twd2h2(1:3,4);
p_h3_wd = Twd2h3(1:3,4);
p_h4_wd = Twd2h4(1:3,4);

% body offset up by h
Tcom2h1_up = [eye(3) [L/2 W/2 h]';
            0 0 0 1];
Tcom2h2_up = [eye(3) [L/2 -W/2 h]';
            0 0 0 1];
Tcom2h3_up = [eye(3) [-L/2 W/2 h]';
            0 0 0 1];
Tcom2h4_up = [eye(3) [-L/2 -W/2 h]';
            0 0 0 1];
Twd2h1_up = Twd2com * Tcom2h1_up;
Twd2h2_up = Twd2com * Tcom2h2_up;
Twd2h3_up = Twd2com * Tcom2h3_up;
Twd2h4_up = Twd2com * Tcom2h4_up;

p_h1_up = Twd2h1_up(1:3,4);
p_h2_up = Twd2h2_up(1:3,4);
p_h3_up = Twd2h3_up(1:3,4);
p_h4_up = Twd2h4_up(1:3,4);

chain1 = [p_h1_wd,p_h2_wd,p_h4_wd,p_h3_wd];
chain2 = [p_h1_wd,p_h2_wd,p_h2_up,p_h1_up];
chain3 = [p_h1_wd,p_h3_wd,p_h3_up,p_h1_up];
chain4 = [p_h3_wd,p_h4_wd,p_h4_up,p_h3_up];
chain5 = [p_h4_wd,p_h2_wd,p_h2_up,p_h4_up];
chain6 = [p_h1_up,p_h2_up,p_h4_up,p_h3_up];

%% inverse kinematics

% --- the main line ---
q = zeros(12,1);
chain_leg = zeros(3,4,4);
for i_leg = 1:4
    if i_leg == 1
        p.sign_L = 1;
        p.sign_d = 1;
    elseif i_leg == 2
        p.sign_L = 1;
        p.sign_d = -1;
    elseif i_leg == 3
        p.sign_L = -1;
        p.sign_d = 1;
    elseif i_leg == 4
        p.sign_L = -1;
        p.sign_d = -1;
    end
    
    q_idx = 3*(i_leg - 1) + (1:3); %3*i_leg-2 : 3*i_leg;
    q(q_idx) = fcn_invKin3(Xd,pf34(:,i_leg),p);
    chain_leg(:,:,i_leg) = legKin(Twd2com,q(q_idx),p);
end
% ---------------------
chain_leg1 = chain_leg(:,:,1);
chain_leg2 = chain_leg(:,:,2);
chain_leg3 = chain_leg(:,:,3);
chain_leg4 = chain_leg(:,:,4);

%% plot

% body
f1 = fill3(chain1(1,:),chain1(2,:),chain1(3,:),body_color,...
           chain2(1,:),chain2(2,:),chain2(3,:),body_color,...
           chain3(1,:),chain3(2,:),chain3(3,:),body_color,...
           chain4(1,:),chain4(2,:),chain4(3,:),body_color,...
           chain5(1,:),chain5(2,:),chain5(3,:),body_color,...
           chain6(1,:),chain6(2,:),chain6(3,:),body_color);
alpha(f1,0.2)

% 
% % legs
% plot3(chain_leg1(1,:),chain_leg1(2,:),chain_leg1(3,:),'linewidth',3,'color',leg_color)
% plot3(chain_leg2(1,:),chain_leg2(2,:),chain_leg2(3,:),'linewidth',3,'color',leg_color)
% plot3(chain_leg3(1,:),chain_leg3(2,:),chain_leg3(3,:),'linewidth',3,'color',leg_color)
% plot3(chain_leg4(1,:),chain_leg4(2,:),chain_leg4(3,:),'linewidth',3,'color',leg_color)
% 
% % feet
% plot3(pf34(1,1),pf34(2,1),pf34(3,1),'o','MarkerFaceColor',leg_color,'MarkerEdgeColor',leg_color)
% plot3(pf34(1,2),pf34(2,2),pf34(3,2),'o','MarkerFaceColor',leg_color,'MarkerEdgeColor',leg_color)
% plot3(pf34(1,3),pf34(2,3),pf34(3,3),'o','MarkerFaceColor',leg_color,'MarkerEdgeColor',leg_color)
% plot3(pf34(1,4),pf34(2,4),pf34(3,4),'o','MarkerFaceColor',leg_color,'MarkerEdgeColor',leg_color)
% 
% 
% GRF
scale = 1e-2;
for i_leg = 1:4
    chain_f = [pf34(:,i_leg),pf34(:,i_leg) + scale * f34(:,i_leg)];
    plot3(chain_f(1,:),chain_f(2,:),chain_f(3,:),'g','linewidth',1.5)
end


end


function chain = legKin(Twd2com,q,p)
    L = p.L;
    W = p.W;
    d = p.d;
    l1 = p.l1;
    l2 = p.l2;
    sign_L = p.sign_L;
    sign_d = p.sign_d;
    
    Tcom2h = [rx(q(1)) [sign_L*L/2 sign_d*W/2 0]';
                0 0 0 1];
    Th2s = [ry(q(2)) [0 sign_d*d 0]';
                0 0 0 1];
    Ts2k = [ry(q(3)) [l1 0 0]';
                0 0 0 1];
    Tk2f = [eye(3) [l2 0 0]';
                0 0 0 1];
    Twd2h = Twd2com * Tcom2h;
    Twd2s = Twd2h * Th2s;
    Twd2k = Twd2s * Ts2k;
    Twd2f = Twd2k * Tk2f;

    p_h_wd = Twd2h(1:3,4);
    p_s_wd = Twd2s(1:3,4);
    p_k_wd = Twd2k(1:3,4);
    p_f_wd = Twd2f(1:3,4);

    chain = [p_h_wd p_s_wd p_k_wd p_f_wd];
end











































