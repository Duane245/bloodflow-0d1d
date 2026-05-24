function dx = FunLeftHeart(t,y)
% parameters
RA = 0.01;
RC = 0.0398;
LS = 0.001025;

RS = 0.8738;
CR = 4.00;

V0 = 15.0;
CS = 2.896;
VD = 5.0;

RM = 0.005;
Emax = 3.0;
Ep = 0.06;          %左心室的被动弹性



T_cycle = 0.8;
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
% 计算时变弹性函数E
E = Ea + Ep;

PLV = Ea*(y(1)-VD) + Ep*(y(1)-V0);
if y(4) > PLV
    DM = 1;
else
    DM = 0;
end

if PLV > y(3)
    DA = 1;
else
    DA = 0;
end

% 计算A和b矩阵
A = [-(Ea+Ep)*DM/RM -DA 0 DM/RM;
    (Ea+Ep)*DA/LS -(DA*RA+RC)/LS -DA/LS 0;
    0 1/CS -1/(RS*CS) 1/(RS*CS);
    (Ea+Ep)*DM/(CR*RM) 0 1/(CR*RS) -(1/(CR*RS)+DM/(CR*RM))];
b = [(E*VD+Ep*(V0-VD))*DM/RM;-(E*VD+Ep*(V0-VD))*DA/LS;
    0;-(E*VD+Ep*(V0-VD))*DM/(CR*RM)];

dx = A*[y(1);y(2);y(3);y(4)]+b;
end
