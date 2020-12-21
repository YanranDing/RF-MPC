function [u_ext,p_ext] = fcn_get_disturbance(t,p)

bz_w = [0 0.5 1 1 0.5 0];

if (t >= 0.5) && (t <= 1.3) % small disturbances
    s_w = (t - 0.5) / 0.8;
    w = polyval_bz(8*bz_w,s_w);
elseif (2.3 <= t) && (t <= 3.1)
    s_w = (t - 2.3) / 0.8;
    w = polyval_bz(22*bz_w,s_w);
else
    w = 0;
end

u_ext = [0;w;0];

p_ext = [p.L/2;p.W/2;p.d];  % external force point in body frame
    
end