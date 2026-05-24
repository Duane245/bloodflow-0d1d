clear
%%
% parameters

CR = 4.4;
CA = 0.08;

RS = 1;
LS = 0.0005;
CS = 1.33;
RC = 0.0398;

R = 0;

RA = 0.001;
RM = 0.005;

Emax = 2.0-0.06;
Ep = 0.06;          %左心室的被动弹性
HR = 75;
%%
step = 0.5e-5;

T_cycle = 60/HR;
t = 0:step:T_cycle
tn = t/(0.2+0.1555*T_cycle);

% 计算主动弹性函数Ea的归一化函数
tn1 = tn/0.7;
tn2 = tn/1.173474;
En1 = tn1.^1.9./(1+tn1.^1.9);
En2 = 1./(1+tn2.^21.9);
En3 = En1.*En2;
En = 1.553174*En3;

clear tn1 tn2 En1 En2 En3
plot(tn,En)

% 计算主动弹性函数Ea
Ea = Emax*En;
plot(t,Ea)
% 计算时变弹性函数E
E = Ea + Ep;
plot(t,E)
%%
C = 1./E;
plot(t,C)
dC = diff(C)./step;
dC = [dC 0];
plot(t,dC)
%%
dC_div_C = -dC./C;
plot(t,dC_div_C)
%%
LVP = 6.9;
LAP = 9.6;
AP = 67;
AOP = 80;
Q = 0;

Y0 = [LVP;LAP;AP;AOP;Q];
%%
% Ea1 = [Ea Ea Ea Ea Ea Ea];
% Ea = Ea1;
% E1 = [E E E E E E];
% E = E1;
C = [C C C C];
dC = [dC dC dC dC];
NT = 4;
%%
time_step = NT*size(t,2);
% time_step = size(t,2)
Y = zeros(time_step,5);
Y(1,:) = Y0';
%%

for n = 1:time_step-1
    
    % 计算A和b矩阵
    A = [-dC(n)/(C(n)+R) 0 0 0 0;
        0 -1/(RS*CR) 1/(RS*CR) 0 0;
        0 1/(RS*CS) -1/(RS*CS) 0 1/CS;
        0 0 0 0 -1/CA;
        0 0 -1/LS 1/LS -RC/LS];
    BMtr = [1/(C(n)+R) -1/(C(n)+R);
        -1/CR 0;
        0 0
        0 1/CA;
        0 0];
    CMtr = [r(Y(n,2)-Y(n,1))/RM;
        r(Y(n,1)-Y(n,4))/RA];
    b = BMtr*CMtr;
    
    [T,Y1]=ode45(@(t,y) FunLeftHeart1(t,y,A,b),...
        [(n-1)*step n*step],Y(n,:)');
    Y(n+1,:) = Y1(end,:);
    
end
%%
% plot(t,Y(:,1))
% plot(t,Y(:,2))
% plot(t,Y(:,3))
% plot(t,Y(:,4))
% plot(t,Y(:,5))

plot(Y(:,1))
plot(Y(:,2))
plot(Y(:,3))
plot(Y(:,4))
plot(Y(:,5))
%%
plot(Y(:,1),'linewidth',2)
hold on
plot(Y(:,2),'linewidth',2)
plot(Y(:,3),'linewidth',2)
hold off
legend('LVP','LAP','AOP')
%%
V0 = 5;
LVV = Y(:,1)./E + V0;
plot(PVV,'linewidth',2)
