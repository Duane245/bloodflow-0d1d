function E = HeartAtrium(t,T0,EA,EB,T_acp,T_arp,t_ac,t_ar)
e_at = zeros(size(t));
idx = find(t<=t_ar +T_arp - T0);
e_at(idx) = 0.5*(1+cos(pi*(t(idx)+T0-t_ar)/T_arp));
idx = intersect(find(t>t_ac),find(t<=t_ac+T_acp));
e_at(idx) = 0.5*(1-cos(pi*(t(idx)-t_ac)/T_acp));
idx = intersect(find(t>t_ac+T_acp),find(t<=T0));
e_at(idx) = 0.5*(1+cos(pi*(t(idx)-t_ar)/T_arp));
E = EA*e_at+EB;
end
