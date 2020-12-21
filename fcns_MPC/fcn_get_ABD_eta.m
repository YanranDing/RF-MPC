function [A,B,D] = fcn_get_ABD_eta(Xt,Ut,p)
% linear dynamics for rotation
% evolution variable is eta

%% parameters
dt = p.Tmpc;

%% unpack
xop = reshape(Xt(1:3),[3,1]);
vop = reshape(Xt(4:6),[3,1]);
Rop = reshape(Xt(7:15),[3,3]);
wop = reshape(Xt(16:18),[3,1]);
pf34 = reshape(Xt(19:30),[3,4]);

%% constants for linear matrices
% [x,v,eta,w,constant]
[Cx_x,Cx_v,Cv_v,Cv_u,Cv_c] = eta_co_xv(Ut,dt,p.mass,p.g);
[CE_eta, CE_w, CE_c] = eta_co_R(Rop,wop,dt);
[Cw_x,Cw_eta,Cw_w, Cw_u, Cw_c] = eta_co_w(xop,Rop,wop,Ut,dt,p.J,pf34);


%% Assemble matrices
    A = [Cx_x, Cx_v, zeros(3,6);
         zeros(3), Cv_v, zeros(3,6);
         zeros(3,6),CE_eta,CE_w;
         Cw_x,zeros(3),Cw_eta,Cw_w];
    B = [zeros(3,12);
         Cv_u;
         zeros(3,12);
         Cw_u];
    D = [zeros(3,1);
         Cv_c;
         CE_c;
         Cw_c];

end

%% Core fcns for constant matrix
function [Cx_x,Cx_v,Cv_v,Cv_u,Cv_c] = eta_co_xv(fop,dt,mass,g)

Cx_x = eye(3);
Cx_v = eye(3) * dt;

Cv_v = eye(3);
Cv_u = dt/mass * [eye(3),eye(3),eye(3),eye(3)];
Cv_c = Cv_u * fop + [0;0;-g] * dt;

end

function [CE_eta, CE_w, CE_c] = eta_co_R(Rop,wop,dt)
% the input arguments are composed of variables at the operating point 
% and parameters

N = fcn_get_N;

%% debugged code
invN = pinv(N);

C_eta = kron(eye(3),Rop*hatMap(wop))*N + kron(eye(3),Rop)*fcn_get_D(wop);
C_w = kron(eye(3),Rop) * N;
C_c = vec(Rop*hatMap(wop)) - kron(eye(3),Rop)*N*wop;

CE_eta = eye(3) + invN * dt * kron(eye(3),Rop') * C_eta;
CE_w = invN * dt * kron(eye(3),Rop') * C_w;
CE_c = invN * dt * kron(eye(3),Rop') * C_c;

end

function [Cw_x,Cw_eta,Cw_w, Cw_u, Cw_c] = eta_co_w(xop,Rop,wop,fop,dt,J,pf)
% the input arguments are composed of variables at the operating point 
% and parameters

N = fcn_get_N;
r1 = pf(:,1) - xop;
r2 = pf(:,2) - xop;
r3 = pf(:,3) - xop;
r4 = pf(:,4) - xop;
Mop = [hatMap(r1) hatMap(r2) hatMap(r3) hatMap(r4)] * fop;

temp_J_w = hatMap(J*wop) - hatMap(wop) * J;
sum_fop = [eye(3),eye(3),eye(3),eye(3)] * fop;

Cx = Rop' * hatMap(sum_fop);
Ceta = fcn_get_F(Mop) * N - temp_J_w * hatMap(wop);
Cw = temp_J_w;
Cu = Rop' * [hatMap(r1),hatMap(r2),hatMap(r3),hatMap(r4)];
Cc = -hatMap(wop)*J*wop + Rop'*Mop - temp_J_w * wop - Cx*xop;

Cw_x = dt*(J\Cx);
Cw_eta = dt*(J\Ceta);
Cw_w = dt*(J\Cw) + eye(3);
Cw_u = dt*(J\Cu);
Cw_c = dt*(J\Cc);

end

%% Aux fcns
function F = fcn_get_F(k)

F = [k', zeros(1,3),zeros(1,3);...
     zeros(1,3),k',zeros(1,3);...
     zeros(1,3),zeros(1,3),k'];
end

function N = fcn_get_N

N = [0 0 0;...
     0 0 1;...
     0 -1 0;...
     0 0 -1;...
     0 0 0;...
     1 0 0;...
     0 1 0;...
     -1 0 0;...
     0 0 0];
end

function D = fcn_get_D(in)

d = in(1);
e = in(2);
f = in(3);
D = [0 0 0;
   e -d 0;
   f 0 -d;
   -e d 0;
   0 0 0;
   0 f -e;
   -f 0 d;
   0 -f e;
   0 0 0];
end




