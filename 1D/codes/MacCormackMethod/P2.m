clear;

L = 350;
N = 1600;
Ct = 1;

dt = 1e-4;
number_of_time_step = 6000;

x = linspace(0,L,N);
dx = x(2)-x(1);

beta = 0.018734e6;
rho = 1.050e3;
A0 = 3.2168;
c0 = sqrt(beta/2/rho*sqrt(A0))*100;

% Q = sin(((x<80).*x)/80*pi);
Q = zeros(1,N);
A = A0*ones(1,N);

W1 = Q./A + 4*ones(1,N)*c0;
W2 = Q./A - 4*ones(1,N)*c0;
% plot(x,W1)
% plot(x,W2)

lamda1 = Q./A + sqrt(beta/2/rho*sqrt(A))*100;
lamda2 = Q./A - sqrt(beta/2/rho*sqrt(A))*100;

W20 = W2(end);
W10 = W1(end);

Rt = 0;

% Cf = 0;
Cf = (2*A0*c0)/2000;
% Cf = 2000/(2*A0*c0);
Cv = 0;
% Cv = 6275;

for time_step = 1:number_of_time_step
    
    t = time_step*dt;
     
    W2n1 = interp1(x,W2,-lamda2(1)*dt);
    W1(1) = -W2n1 + 2*InflowFuncation(t)/A(1);
%     W1(1) = interp1(x,W1, lamda1(1)*dt);
    W2(1) = W2n1;
    
    W1n1 = interp1(x,W1,L-lamda1(end)*dt);
    W1(N) = W1n1;
    W2(N) = W20 - Rt*(W1n1-W10);
    
    W1i = W1;
    W2i = W2;
    
    [W1,W2] = PreStep(W1,W2,lamda1,lamda2,Q,A,Cf,Cv,dx,dt,N,A0);
    
    [A,Q] = WtoAQ(W1,W2,rho,beta);
    
    [W1,W2] = CorStep(W1i,W2i,W1,W2,lamda1,lamda2,Q,A,Cf,Cv,dx,dt,N,A0);
    
    [A,Q] = WtoAQ(W1,W2,rho,beta);
    
    lamda1 = Q./A + sqrt(beta/2/rho*sqrt(A))*100;
    lamda2 = Q./A - sqrt(beta/2/rho*sqrt(A))*100;
end
% plot(x,W1)
% plot(x,W2)
% plot(x,A)
plot(x,Q)
% plot(x,lamda1)
% plot(x,lamda2)
%%
Cf = (2*A0*c0)/2000;
% Cf = 2000/(2*A0*c0);
DashLine = exp(-Cf*x/(2*A0*c0));
plot(x,DashLine)
hold on 
plot(x,Q)
hold off
