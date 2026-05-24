%% 0D boundary

clear;

L = 40;
N = 4001;
Ct = 1;

dt = 1e-5;
number_of_time_step = 50000;

x = linspace(0,L,N);
dx = x(2)-x(1);

W1 = zeros(1,N);
W2 = zeros(1,N);


% beta = 0.018734e6;
h = 0.15;
E = 0.4e6;
A0 = 3.14;
beta = 4/3*sqrt(pi)*h*E/A0;
rho = 1.050e3;

c0 = sqrt(beta/2/rho*sqrt(A0))*100;

W1 = 4*ones(1,N)*c0;
W2 = -4*ones(1,N)*c0;

A = (W1-W2).^4/1024*(rho/beta)^2*1e-8;
Q = A.*(W1+W2)/2;

lamda1 = Q./A + sqrt(beta/2/rho*sqrt(A))*100;
lamda2 = Q./A - sqrt(beta/2/rho*sqrt(A))*100;

W20 = W2(end);
W10 = W1(end);

Rt = 0.8;

Cf = 0;
% Cf = (2*A0*c0)/2000;
% Cf = 2000/(2*A0*c0);
Cv = 0;
% Cv = 6275;

% R = 0;
% C = 6.3;
% L = 1;
C = 6.3e-3;
L_inductance = 1e-2;    %
pout = 0;

pc = 0;

Z0 = rho*c0/A0;

% R = 1e4;
% R = Z0 *1e-4;
R = 189;

for time_step = 1:number_of_time_step
    
    t = time_step*dt;
    
    %%%%%%%%%%%%%%%%%%
    UL = Q(end)/A(end);
    AL = A(end);
    F = @(X) ( R*( UL + 4*sqrt(beta/2/rho)*(AL^(1/4))*100) * X - 4*R*sqrt(beta/2/rho)*(X^(1/4))*100*X - beta*(sqrt(X)-sqrt(A0)) + pout);
    X0 = A(end);
    options = optimoptions('fsolve','Display','iter');
    [X,fval] = fsolve(F,X0,options);
    
    Xcup(time_step) = X;
    
    Axing = X;
    PAxing = beta*(sqrt(Axing)-sqrt(A0));
    Uxing = (PAxing - pout)/Axing/R;
    
    UR = 2*Uxing-UL;
    Q(end) = AL*UR;
    %%%%%%%%%%%%%%%%%%%%
    
    W2n1 = interp1(x,W2,-lamda2(1)*dt);
    W1(1) = -W2n1 + 2*InflowFuncation(t)/A(1);
    W2(1) = W2n1;
    
    W1n1 = interp1(x,W1,L-lamda1(end)*dt);
    W1(N) = W1n1;
%     W2(N) = W20 - Rt*(W1n1-W10);
%     W1(N) = Q(end)/A(end) + 4*sqrt(beta/2/rho*sqrt(A(end)))*100;
    W2(N) = Q(end)/A(end) - 4*sqrt(beta/2/rho*sqrt(A(end)))*100;
    
    W1i = W1;
    W2i = W2;
    
    W1(2:N-1) = W1(2:N-1) - (dt/dx)*(W1(3:N)-W1(2:N-1)).*lamda1(2:N-1);
    W2(2:N-1) = W2(2:N-1) - (dt/dx)*(W2(3:N)-W2(2:N-1)).*lamda2(2:N-1) - dt*(Cf/A0)*Q(2:N-1);
    
    A = (W1-W2).^4/1024*(rho/beta)^2*1e-8;
    Q = A.*(W1+W2)/2;
    
    W1(2:N-1) = (W1(2:N-1)+W1i(2:N-1))/2 - (dt/dx/2)*(W1(2:N-1)-W1(1:N-2)).*lamda1(2:N-1);
    W2(2:N-1) = (W2(2:N-1)+W2i(2:N-1))/2 - (dt/dx/2)*(W2(2:N-1)-W2(1:N-2)).*lamda2(2:N-1) - (dt/2)*(Cf/A0)*Q(2:N-1);
    
    A = (W1-W2).^4/1024*(rho/beta)^2*1e-8;
    Q = A.*(W1+W2)/2;
    
    lamda1 = Q./A + sqrt(beta/2/rho*sqrt(A))*100;
    lamda2 = Q./A - sqrt(beta/2/rho*sqrt(A))*100;
    
    Qcapture1(time_step) = Q(floor(N/2));
    
    Acupture1(time_step) = A(floor(N/2));
    
    Ucapture1(time_step) = Q(floor(N/2))/A(floor(N/2));
end
plot(x,W1)
plot(x,W2)
% plot(x,A)
plot(x,Q)
% plot(x,lamda1)
% plot(x,lamda2)

%%
% plot((1:number_of_time_step)*dt,Qcapture1)
% plot((1:number_of_time_step)*dt,Ucapture1/max(Ucapture1))

P = beta*(sqrt(Acupture1)-sqrt(A0));
%%
plot((1:number_of_time_step)*dt,P/max(P),"LineWidth",1.5)
xlabel("t(s)")
ylabel("P/Pmax")
ylim([-1,1])
title("single resistance model (R = 189)")
% title("CR model")
