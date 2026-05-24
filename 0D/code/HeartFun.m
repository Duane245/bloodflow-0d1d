function [t,e_lv] = HeartFun(E_lva,E_lvb,tau_1,tau_2,t_ee,t_r,step_size)
t = 0:step_size:t_r;
idx = find(t<t_ee);
e_lv = zeros(size(t));
e_lv(idx) = E_lva*(1-exp(-t(idx)/tau_1))+E_lvb;
idx = find(t>=t_ee);
e_lv(idx) = (E_lva*(1-exp(-t(idx(1))/tau_1))+E_lvb ...
    -E_lvb)*exp(-(t(idx)-t_ee)/tau_2) + E_lvb;
end
