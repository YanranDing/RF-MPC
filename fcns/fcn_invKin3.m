function q = fcn_invKin3(X,pf,p)
pf = reshape(pf,length(pf),1);
% states
pcom = X(1:3);
R = reshape(X(7:15),[3,3]);

% parameters
L = p.L;
W = p.W;
sign_L = p.sign_L;
sign_d = p.sign_d;
    
%%
Twd2com = [R pcom;
           0 0 0 1];
Tcom2h = [eye(3) [sign_L*L/2 sign_d*W/2 0]';
            0 0 0 1];
Twd2h = Twd2com * Tcom2h;
p_f_wd = pf(1:3);
p_h_wd = Twd2h(1:3,4);
        
% h1 to f1 in world frame
p_h2f_wd = p_f_wd - p_h_wd;

% h1 to f1 in body frame
p_h2f_b = R' * p_h2f_wd;

%% inverse kinematics
q = invKin(p_h2f_b,p);

end


function q = invKin(p_h2f_b,p)
l1 = p.l1;
l2 = p.l2;
d = p.d;
sign_d = p.sign_d;

if size(p_h2f_b,1) == 1
    p_h2f_b = p_h2f_b';
end
vp = p_h2f_b;

%% q1 
vpyz = vp(2:3);
ryz = norm(vpyz,2);
a = asin(vp(2)/ryz);
b = asin(d/ryz);
q1 = a - sign_d * b;

%% q2
r = norm(vp,2);
vd = sign_d * [0 d*cos(q1) d*sin(q1)]';
rf = norm(vp - vd,2);
vf = rx(q1)' * (vp - vd);
if vf(3) <= 0
    a = acos(vf(1)/rf);
else
    if vf(1) >=0
        a = -asin(vf(3)/rf);
    else
        a = pi + asin(vf(3)/rf);
    end
end
cosb = (l1^2 + rf^2 - l2^2)/(2*l1*rf);
b = acos(cosb);
q2 = a + b;

%% q3
cosc = (l1^2+l2^2-rf^2)/(2*l1*l2);
c = acos(cosc);
q3 = -(pi - c);

q = [q1 q2 q3];

end



