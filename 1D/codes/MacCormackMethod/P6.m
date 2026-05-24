clear;

number_of_vessel = 55;

initial_L;
% N = 1600;
% Ct = 1;

dt = 4e-7;
number_of_time_step = 1e6;

% x = linspace(0,L(1),N);
% dx = x(2)-x(1);

initial_beta;
rho = 1.050e3;

initial_A0;

for i = 1:number_of_vessel
    c0_str = ['c0(' num2str(i) ...
        ')=sqrt((beta(' num2str(i) ')/2/rho)*sqrt(A0(' num2str(i) ')))*100;'];
    eval(c0_str);
end


for i = 1:number_of_vessel
    Rl(i) = L(i)/c0(i);
end
dx = 0.1;
for i = 1:number_of_vessel
    N(i) = round(L(i)/dx + 1);
end

for i = 1:number_of_vessel
    x{i} = 0:dx:L(i);
end
%%

for i = 1:number_of_vessel
    W1_str = ['W1.artery' num2str(i) '=4*ones(1,' num2str(N(i)) ')*c0(' num2str(i) ');'];
    eval(W1_str);
end

for i = 1:number_of_vessel
    W2_str = ['W2.artery' num2str(i) '=-4*ones(1,' num2str(N(i)) ')*c0(' num2str(i) ');'];
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

%%
for i = 1:number_of_vessel
    str1 = ['W10(' num2str(i) ') = W1.artery' num2str(i) '(end);'];
    eval(str1);
    str2 = ['W20(' num2str(i) ') = W2.artery' num2str(i) '(end);'];
    eval(str2);
end
%%
initial_Rt;

Cf = 0;
% Cf = (2*A0*c0)/2000;
% Cf = 2000/(2*A0*c0);

initial_Cv;
% Cv = 6275;
%%
idxOutflowArtery = [6 8 10 11 12 13 16 17 20 22 24 25 26 31 32 33 34 36 38 40 45 47 48 49 51 53 54 55];
CB = [1 2 3;
    3 4 5;
    5 12 13;
    4 6 7;
    7 8 9;
    9 10 11;
    2 14 15;
    15 16 17;
    14 18 19;
    19 20 21;
    21 22 23;
    23 24 25;
    18 26 27;
    27 28 29;
    29 30 31;
    30 32 33;
    28 34 35;
    35 36 37;
    37 38 39;
    39 40 41;
    41 42 43;
    43 50 51;
    50 52 53;
    52 54 55;
    42 44 45;
    44 46 47;
    46 48 49];
number_of_CB = size(CB,1);
%%

for time_step = 1:number_of_time_step
    t = time_step*dt;
    
    % Parent inflow boundary condition (n time step W2 computing n+1 timestep W1)
    W2n1 = interp1(x{1},W2.artery1,-lamda2.artery1(1)*dt);
    W1.artery1(1) = -W2n1 + 2*500*InflowFuncation(t)/A.artery1(1);
    W2.artery1(1) = W2n1;
    
    % outflow boundary condition
    for i = idxOutflowArtery
        str1 = ['W1n1 = interp1(x{' num2str(i) '},W1.artery' num2str(i) ',L(' num2str(i) ')-lamda1.artery' num2str(i) '(end)*dt);'];
        eval(str1);
        str2 = ['W1.artery' num2str(i) '(N(' num2str(i) ')) = W1n1;'];
        eval(str2);
        str3 = ['W2.artery' num2str(i) '(N(' num2str(i) ')) = W20(' num2str(i) ') - Rt(' num2str(i) ')*(W1n1-W10(' num2str(i) '));'];
        eval(str3);
    end
    
    for i = 1:number_of_CB
        ConnectionBoudary = CB(i,:);
        BoundaryCondition;
    end
    
    W1i = W1;
    W2i = W2;
    
    for i = 1:number_of_vessel
        str1 = ['[W1.artery' num2str(i) ',W2.artery' num2str(i) ']=PreStep(W1.artery' num2str(i) ',W2.artery' num2str(i) ',lamda1.artery' num2str(i) ',lamda2.artery' num2str(i) ',Q.artery' num2str(i) ',A.artery' num2str(i) ',Cf,Cv(' num2str(i) '),dx,dt,N(' num2str(i) '),A0(' num2str(i) '));'];
        eval(str1);
        str2 = ['[A.artery' num2str(i) ',Q.artery' num2str(i) '] = WtoAQ(W1.artery' num2str(i) ',W2.artery' num2str(i) ',rho,beta(' num2str(i) '));';];
        eval(str2);
        str3 = ['[W1.artery' num2str(i) ',W2.artery' num2str(i) '] = CorStep(W1i.artery' num2str(i) ',W2i.artery' num2str(i) ',W1.artery' num2str(i) ',W2.artery' num2str(i) ',lamda1.artery' num2str(i) ',lamda2.artery' num2str(i) ',Q.artery' num2str(i) ',A.artery' num2str(i) ',Cf,Cv(' num2str(i) '),dx,dt,N(' num2str(i) '),A0(' num2str(i) '));'];
        eval(str3);
        str4 = ['[A.artery' num2str(i) ',Q.artery' num2str(i) '] = WtoAQ(W1.artery' num2str(i) ',W2.artery' num2str(i) ',rho,beta(' num2str(i) '));'];
        eval(str4);
        str5 = ['lamda1.artery' num2str(i) '= Q.artery' num2str(i) './A.artery' num2str(i) '+ sqrt(beta(' num2str(i) ')/2/rho*sqrt(A0(' num2str(i) ')))*100;'];
        eval(str5);
        str6 = ['lamda2.artery' num2str(i) '= Q.artery' num2str(i) './A.artery' num2str(i) '- sqrt(beta(' num2str(i) ')/2/rho*sqrt(A0(' num2str(i) ')))*100;'];
        eval(str6);
    end
    
    captureQ1(time_step) = Q.artery1(20);
    captureQ2(time_step) = Q.artery8(118);
    captureQ3(time_step) = Q.artery37(6);
    captureQ4(time_step) = Q.artery54(162);
    
    captureA1(time_step) = A.artery1(20);
    captureA2(time_step) = A.artery8(118);
    captureA3(time_step) = A.artery37(6);
    captureA4(time_step) = A.artery54(162);
    
%     captureA1(time_step) = A.artery1(120);
%     captureA2(time_step) = A.artery2(100);
end
%%
% plot(x{1},W1.artery1)
% plot(x{2},W1.artery2)
% plot(x{3},W1.artery3)
% plot(x{1},W2.artery1)
% plot(x{2},W2.artery2)
% plot(x{3},W2.artery3)
% plot(x,A)
% plot(x{1},Q.artery1)
% plot(x{2},Q.artery2)
% plot(x{3},Q.artery3)
% plot(x,lamda1)
% plot(x,lamda2)
%%
% PA = beta(1)*((sqrt(captureA1))-(sqrt(A0(1))));
% PB = beta(2)*((sqrt(captureA2))-(sqrt(A0(2))));
% figure
% plot(PA)
% hold on
% plot(PB)
% hold off
%%
% PNormal = max(PA);
% PA = beta(1)/PNormal*((sqrt(captureA1))-(sqrt(A0(1))));
% PB = beta(2)/PNormal*((sqrt(captureA2))-(sqrt(A0(2))));
% figure
% plot(PA)
% hold on
% plot(PB)
% hold off
