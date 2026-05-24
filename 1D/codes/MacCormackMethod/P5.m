clear;

number_of_vessel = 3;

L = 250;
N = 1600;
Ct = 1;

dt = 1e-4;
number_of_time_step = 12000;

x = linspace(0,L,N);
dx = x(2)-x(1);

beta(1) = 0.023633e6;
beta(2) = 0.063021e6;
beta(3) = 0.063021e6;
rho = 1.050e3;

A0(1) = 4;
A0(2) = 1.5;
A0(3) = 1.5;

for i = 1:number_of_vessel
    c0_str = ['c0(' num2str(i) ...
        ')=sqrt((beta(' num2str(i) ')/2/rho)*sqrt(A0(' num2str(i) ')))*100;'];
    eval(c0_str);
end

for i = 1:number_of_vessel
    W1_str = ['W1.artery' num2str(i) '=4*ones(1,' num2str(N) ')*c0(' num2str(i) ');'];
    eval(W1_str);
end

for i = 1:number_of_vessel
    W2_str = ['W2.artery' num2str(i) '=-4*ones(1,' num2str(N) ')*c0(' num2str(i) ');'];
    eval(W2_str);
end

for i = 1:number_of_vessel
    A_str = ['A.artery' num2str(i) '=(W1.artery' num2str(i) '-W2.artery' num2str(i) ').^4/1024*(rho/beta(' num2str(i) '))^2*1e-8;'];
    eval(A_str);
end

for i = 1:number_of_vessel
    Q_str = ['Q.artery' num2str(i) '=A.artery' num2str(i) '.*(W1.artery' num2str(i) '+W2.artery' num2str(i) ')/2;'];
    eval(Q_str);
end

for i = 1:number_of_vessel
    lamda1_str = ['lamda1.artery' num2str(i) '=Q.artery' num2str(i) './A.artery' num2str(i) '+sqrt(beta(' num2str(i) ')/2/rho*sqrt(A.artery' num2str(i) '))*100;'];
    eval(lamda1_str);
end

for i = 1:number_of_vessel
    lamda2_str = ['lamda2.artery' num2str(i) '=Q.artery' num2str(i) './A.artery' num2str(i) '-sqrt(beta(' num2str(i) ')/2/rho*sqrt(A.artery' num2str(i) '))*100;'];
    eval(lamda2_str);
end

clear W1_str W2_str c0_str A_str Q_str lamda1_str lamda2_str

W1parent0 = W1.artery1(end);
W2parent0 = W2.artery1(end);
W1daughter10 = W1.artery2(end);
W2daughter10 = W2.artery2(end);
W1daughter20 = W1.artery3(end);
W2daughter20 = W2.artery3(end);
%%
Rt = 0;

Cf = 0;
% Cf = (2*A0*c0)/2000;
% Cf = 2000/(2*A0*c0);
Cv = 0;
% Cv = 6275;
%%

for time_step = 1:number_of_time_step
    t = time_step*dt;
    
    % Parent inflow boundary condition (n time step W2 computing n+1 timestep W1)
    W2n1 = interp1(x,W2.artery1,-lamda2.artery1(1)*dt);
    W1.artery1(1) = -W2n1 + 2*InflowFuncation(t)/A.artery1(1);
    W2.artery1(1) = W2n1;
    
    % daughter1 outflow boundary condition
    W1daughter1n1 = interp1(x,W1.artery2,L-lamda1.artery2(end)*dt);
    W1.artery2(N) = W1daughter1n1;
    W2.artery2(N) = W2daughter10 - Rt*(W1daughter1n1-W1daughter10);
    
    % daughter2 outflow boundary condition
    W1daughter2n1 = interp1(x,W1.artery3,L-lamda1.artery3(end)*dt);
    W1.artery3(N) = W1daughter2n1;
    W2.artery3(N) = W2daughter20 - Rt*(W1daughter2n1-W1daughter20);
    
    W1parentn1 = interp1(x,W1.artery1,L-lamda1.artery1(end)*dt);
    W2daughter1n1 = interp1(x,W2.artery2,-lamda2.artery2(1)*dt);
    W2daughter2n1 = interp1(x,W2.artery3,-lamda2.artery3(1)*dt);
    
    F = @(X) [abs(X(2))-abs(X(4))-abs(X(6));
        0.5*rho*(X(2)/X(1))^2+abs(beta(1)*((sqrt(abs(X(1))))-(sqrt(A0(1))))*1e4)-0.5*rho*(X(4)/X(3))^2-abs(beta(2)*((sqrt(abs(X(3))))-(sqrt(A0(2))))*1e4);
        0.5*rho*(X(2)/X(1))^2+abs(beta(1)*((sqrt(abs(X(1))))-(sqrt(A0(1))))*1e4)-0.5*rho*(X(6)/X(5))^2-abs(beta(3)*((sqrt(abs(X(5))))-(sqrt(A0(3))))*1e4);
        -abs(X(2)/X(1))-4*sqrt(beta(1)/2/rho*sqrt(abs(X(1))))*100+W1parentn1;
        -abs(X(4)/X(3))+4*sqrt(beta(2)/2/rho*sqrt(abs(X(3))))*100+W2daughter1n1;
        -abs(X(6)/X(5))+4*sqrt(beta(3)/2/rho*sqrt(abs(X(5))))*100+W2daughter2n1];
    X0 = [A.artery1(end);Q.artery1(end);A.artery2(1);Q.artery2(1);A.artery3(1);Q.artery3(1)];
    options = optimoptions('fsolve','Display','iter');
    [X,fval] = fsolve(F,X0,options);
    
    A.artery1(end) = abs(X(1));
    Q.artery1(end) = abs(X(2));
    A.artery2(1) = abs(X(3));
    Q.artery2(1) = abs(X(4));
    A.artery3(1) = abs(X(5));
    Q.artery3(1) = abs(X(6));
    
    %     X_temp4(time_step) = X(4);
    %     X_temp6(time_step) = X(6);
    %     X_temp2(time_step) = X(2);
    %
    %     X_temp1(time_step) = X(1);
    %     X_temp3(time_step) = X(3);
    %     X_temp5(time_step) = X(5);
    
    % Parent outflow boundary condition
    W1.artery1(N) = Q.artery1(end)./A.artery1(end) + 4*sqrt(beta(1)/2/rho*sqrt(abs(A.artery1(end))))*100;
    W2.artery1(N) = Q.artery1(end)./A.artery1(end) - 4*sqrt(beta(1)/2/rho*sqrt(abs(A.artery1(end))))*100;
    
    % daughter1 inflow boundary condition (n time step W2 computing n+1 timestep W1)
    W1.artery2(1) = Q.artery2(1)./A.artery2(1) + 4*sqrt(beta(2)/2/rho*sqrt(abs(A.artery2(1))))*100;
    W2.artery2(1) = Q.artery2(1)./A.artery2(1) - 4*sqrt(beta(2)/2/rho*sqrt(abs(A.artery2(1))))*100;
    
    % daughter2 inflow boundary condition (n time step W2 computing n+1 timestep W1)
    W1.artery3(1) = Q.artery3(1)./A.artery3(1) + 4*sqrt(beta(3)/2/rho*sqrt(abs(A.artery3(1))))*100;
    W2.artery3(1) = Q.artery3(1)./A.artery3(1) - 4*sqrt(beta(3)/2/rho*sqrt(abs(A.artery3(1))))*100;
    
    
    W1i = W1;
    W2i = W2;
    
    % updating 1 vessel
    [W1.artery1,W2.artery1] = PreStep(W1.artery1,W2.artery1,lamda1.artery1,lamda2.artery1,Q.artery1,A.artery1,Cf,Cv,dx,dt,N,A0(1));
    
    [A.artery1,Q.artery1] = WtoAQ(W1.artery1,W2.artery1,rho,beta(1));
    
    [W1.artery1,W2.artery1] = CorStep(W1i.artery1,W2i.artery1,W1.artery1,W2.artery1,lamda1.artery1,lamda2.artery1,Q.artery1,A.artery1,Cf,Cv,dx,dt,N,A0(1));
    
    [A.artery1,Q.artery1] = WtoAQ(W1.artery1,W2.artery1,rho,beta(1));
    
    lamda1.artery1 = Q.artery1./A.artery1 + sqrt(beta(1)/2/rho*sqrt(A0(1)))*100;
    lamda2.artery1 = Q.artery1./A.artery1 - sqrt(beta(1)/2/rho*sqrt(A0(1)))*100;
    
    % updating 2 vessel
    [W1.artery2,W2.artery2] = PreStep(W1.artery2,W2.artery2,lamda1.artery2,lamda2.artery2,Q.artery2,A.artery2,Cf,Cv,dx,dt,N,A0(2));
    
    [A.artery2,Q.artery2] = WtoAQ(W1.artery2,W2.artery2,rho,beta(2));
    
    [W1.artery2,W2.artery2] = CorStep(W1i.artery2,W2i.artery2,W1.artery2,W2.artery2,lamda1.artery2,lamda2.artery2,Q.artery2,A.artery2,Cf,Cv,dx,dt,N,A0(2));
    
    [A.artery2,Q.artery2] = WtoAQ(W1.artery2,W2.artery2,rho,beta(2));
    
    lamda1.artery2 = Q.artery2./A.artery2 + sqrt(beta(2)/2/rho*sqrt(A0(2)))*100;
    lamda2.artery2 = Q.artery2./A.artery2 - sqrt(beta(2)/2/rho*sqrt(A0(2)))*100;
    
    % updating 3 vessel
    [W1.artery3,W2.artery3] = PreStep(W1.artery3,W2.artery3,lamda1.artery3,lamda2.artery3,Q.artery3,A.artery3,Cf,Cv,dx,dt,N,A0(3));
    
    [A.artery3,Q.artery3] = WtoAQ(W1.artery3,W2.artery3,rho,beta(3));
    
    [W1.artery3,W2.artery3] = CorStep(W1i.artery3,W2i.artery3,W1.artery3,W2.artery3,lamda1.artery3,lamda2.artery3,Q.artery3,A.artery3,Cf,Cv,dx,dt,N,A0(3));
    
    [A.artery3,Q.artery3] = WtoAQ(W1.artery3,W2.artery3,rho,beta(3));
    
    lamda1.artery3 = Q.artery3./A.artery3 + sqrt(beta(3)/2/rho*sqrt(A0(3)))*100;
    lamda2.artery3 = Q.artery3./A.artery3 - sqrt(beta(3)/2/rho*sqrt(A0(3)))*100;
    
    
    
    
    captureA1(time_step) = A.artery1(800);
    captureA2(time_step) = A.artery2(200);
end
%%
str1 = ['[W1.artery' num2str(i) ',W2.artery' num2str(i) ']='...
    'PreStep(W1.artery' num2str(i) ',W2.artery' num2str(i) ...
    ',lamda1.artery' num2str(i) ',lamda2.artery' num2str(i) ...
    ',Q.artery' num2str(i) ',A.artery' num2str(i) ...
    ',Cf,Cv,dx,dt,N,A0(' num2str(i) '));']
% eval(str1); 
str2 = ['[A.artery' num2str(i) ',Q.artery' num2str(i) '] = WtoAQ(W1.artery' num2str(i) ',W2.artery' num2str(i) ',rho,beta(' num2str(i) '));';]
str3 = ['[W1.artery' num2str(i) ',W2.artery' num2str(i) '] = CorStep(W1i.artery' num2str(i) ',W2i.artery' num2str(i) ',W1.artery' num2str(i) ',W2.artery' num2str(i) ',lamda1.artery' num2str(i) ',lamda2.artery' num2str(i) ',Q.artery' num2str(i) ',A.artery' num2str(i) ',Cf,Cv,dx,dt,N,A0(' num2str(i) '));']
str4 = ['[A.artery' num2str(i) ',Q.artery' num2str(i) '] = WtoAQ(W1.artery' num2str(i) ',W2.artery' num2str(i) ',rho,beta(' num2str(i) '));']
str5 = ['lamda1.artery' num2str(i) '= Q.artery' num2str(i) './A.artery' num2str(i) '+ sqrt(beta(' num2str(i) ')/2/rho*sqrt(A0(' num2str(i) ')))*100;']
str6 = ['lamda2.artery' num2str(i) '= Q.artery' num2str(i) './A.artery' num2str(i) '- sqrt(beta(' num2str(i) ')/2/rho*sqrt(A0(' num2str(i) ')))*100;']
%%
plot(x,W1.artery1)
plot(x,W1.artery2)
plot(x,W1.artery3)
plot(x,W2.artery1)
plot(x,W2.artery2)
plot(x,W2.artery3)
% plot(x,A)
plot(x,Q.artery1)
plot(x,Q.artery2)
plot(x,Q.artery3)
% plot(x,lamda1)
% plot(x,lamda2)
%%
PA = beta(1)*((sqrt(captureA1))-(sqrt(A0(1))));
PB = beta(2)*((sqrt(captureA2))-(sqrt(A0(2))));
figure
plot(PA)
hold on
plot(PB)
hold off
%%
PNormal = max(PA);
PA = beta(1)/PNormal*((sqrt(captureA1))-(sqrt(A0(1))));
PB = beta(2)/PNormal*((sqrt(captureA2))-(sqrt(A0(2))));

figure
plot((1:number_of_time_step)*dt,PA,'lineWidth',2)
xlabel("Time(s)")
ylabel("Pressure")
hold on
plot((1:number_of_time_step)*dt,PB,'lineWidth',2)
hold off
