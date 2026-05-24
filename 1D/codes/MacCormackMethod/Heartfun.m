function E = Heartfun(t)
Emax = 2.0-0.06;
Ep = 0.06;          %左心室的被动弹性
HR = 75;
%%

T_cycle = 60/HR;
% t = 0:step:T_cycle;
tn = t/(0.2+0.1555*T_cycle);

% 计算主动弹性函数Ea的归一化函数
tn1 = tn/0.7;
tn2 = tn/1.173474;
En1 = tn1.^1.9./(1+tn1.^1.9);
En2 = 1./(1+tn2.^21.9);
En3 = En1.*En2;
En = 1.553174*En3;

% 计算主动弹性函数Ea
Ea = Emax*En;
% plot(t,Ea)
% 计算时变弹性函数E
E = Ea + Ep;
% plot(t,E)
end
