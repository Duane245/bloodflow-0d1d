clear;

L = 350;
N = 1600;
Ct = 1;

dt = 1e-6;
number_of_time_step = 200000;

x = linspace(0,L,N);
dx = x(2)-x(1);

beta = 0.018734e6;
rho = 1.050e3;
A0 = 3.2168;
c0 = sqrt(beta/2/rho*sqrt(A0))*100;

Q = sin(((x<80).*x)/80*pi);
% Q = zeros(1,N);
plot(x,Q)
A = A0*ones(1,N);

W1 = Q./A + 4*ones(1,N)*c0;
W2 = Q./A - 4*ones(1,N)*c0;
plot(x,W1)
plot(x,W2)

% A = (W1-W2).^4/1024*(rho/beta)^2*1e-8;
% Q = A.*(W1+W2)/2;

lamda1 = Q./A + sqrt(beta/2/rho*sqrt(A))*100;
lamda2 = Q./A - sqrt(beta/2/rho*sqrt(A))*100;

W20 = W2(end);
W10 = W1(end);

Rt = 0;

Cf = 0;
% Cf = (2*A0*c0)/2000;
% Cf = 2000/(2*A0*c0);
% Cv = 0;
Cv = 6275;

for time_step = 1:number_of_time_step
    
    t = time_step*dt;
     
    W2n1 = interp1(x,W2,-lamda2(1)*dt);
%     W1(1) = -W2n1 + 2*InflowFuncation(t)/A(1);
    W1(1) = interp1(x,W1, lamda1(1)*dt);
    W2(1) = W2n1;
    
    W1n1 = interp1(x,W1,L-lamda1(end)*dt);
    W1(N) = W1n1;
    W2(N) = W20 - Rt*(W1n1-W10);
    
    W1i = W1;
    W2i = W2;
    
    dQ2(2:N-1) = (Q(3:N)-2*Q(2:N-1)+Q(1:N-2))/(dx^2);
    
    
    W1(2:N-1) = W1(2:N-1) - (dt/dx)*(W1(3:N)-W1(2:N-1)).*lamda1(2:N-1)+ (dt*Cv)*dQ2(2:N-1)./A(2:N-1);
    W2(2:N-1) = W2(2:N-1) - (dt/dx)*(W2(3:N)-W2(2:N-1)).*lamda2(2:N-1) - dt*(Cf/A0)*Q(2:N-1) + (dt*Cv)*dQ2(2:N-1)./A(2:N-1);
    
    A = (W1-W2).^4/1024*(rho/beta)^2*1e-8;
    Q = A.*(W1+W2)/2;
    
    dQ2(2:N-1) = (Q(3:N)-2*Q(2:N-1)+Q(1:N-2))/(dx^2);
    
    
    W1(2:N-1) = (W1(2:N-1)+W1i(2:N-1))/2 - (dt/dx/2)*(W1(2:N-1)-W1(1:N-2)).*lamda1(2:N-1)+ dt/2*Cv*dQ2(2:N-1)./A(2:N-1);
    W2(2:N-1) = (W2(2:N-1)+W2i(2:N-1))/2 - (dt/dx/2)*(W2(2:N-1)-W2(1:N-2)).*lamda2(2:N-1) - (dt/2)*(Cf/A0)*Q(2:N-1) + dt/2*Cv*dQ2(2:N-1)./A(2:N-1);
    
    A = (W1-W2).^4/1024*(rho/beta)^2*1e-8;
    Q = A.*(W1+W2)/2;
    
    lamda1 = Q./A + sqrt(beta/2/rho*sqrt(A))*100;
    lamda2 = Q./A - sqrt(beta/2/rho*sqrt(A))*100;
end
plot(x,W1)
plot(x,W2)
% plot(x,A)
plot(x,Q)
% plot(x,lamda1)
% plot(x,lamda2)
%%
% Cf = (2*A0*c0)/2000;
% % Cf = 2000/(2*A0*c0);
% DashLine = exp(-Cf*x/(2*A0*c0));
% plot(x,DashLine)
% hold on 
% plot(x,Q)
% hold off
