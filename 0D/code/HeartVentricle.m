function E = HeartVentricle(t,T0,EA,EB,T_vcp,T_vrp,t_vc,t_vr)
e_vt = zeros(size(t));
idx = find(t<=T_vcp);
e_vt(idx) = 0.5*(1-cos(pi*t(idx)/T_vcp));
idx = intersect(find(t>T_vcp),find(t<=T_vcp+T_vrp));
e_vt(idx) = 0.5*(1+cos(pi*(t(idx)-T_vcp)/T_vrp));
E = EA*e_vt+EB;
end
