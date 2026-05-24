%% 使用 MacCormack + Crank-Nicolson 混合方法，计算人体 55 条动脉网络

clear;

number_of_vessel = 55;

L = zeros(number_of_vessel,1);
initial_L;

N = L*10 + 1;
%%
dt = 5e-5;
t = 0;
number_of_time_step = 16000*10;

x = cell(number_of_vessel,1);

for i = 1:number_of_vessel
    x{i} = linspace(0,L(i),N(i));
end

dx = 0.1;
%%
beta = zeros(number_of_vessel,1);
initial_beta;

A0 = zeros(number_of_vessel,1);
initial_A0;

rho = 1.050e3;

c0 = zeros(number_of_vessel,1);
for i = 1:number_of_vessel
    c0(i) = sqrt(beta(i)/2/rho*sqrt(A0(i)))*100;
end

Cf = zeros(number_of_vessel,1);

Cv = zeros(number_of_vessel,1);
initial_Cv;

Rt = zeros(number_of_vessel,1);
initial_Rt;

A = cell(number_of_vessel,1);
Q = cell(number_of_vessel,1);
for i = 1:number_of_vessel
    A{i} = A0(i)*ones(size(x{i}));
    Q{i} = zeros(size(x{i}));
end
%%
[E,a,M,K,MA] = ComputingCrankNicolsonCoffecients(number_of_vessel,dt,dx,Cv,c0,N);

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
captureQ1 = zeros(1,number_of_time_step);
captureQ2 = zeros(1,number_of_time_step);
captureQ3 = zeros(1,number_of_time_step);
captureQ4 = zeros(1,number_of_time_step);

captureA1 = zeros(1,number_of_time_step);
captureA2 = zeros(1,number_of_time_step);
captureA3 = zeros(1,number_of_time_step);
captureA4 = zeros(1,number_of_time_step);
%%
for time_step = 1:number_of_time_step
    t = time_step*dt;
    
    init_A = A;
    init_Q = Q;
    
    for i = 1:number_of_CB
        ConnectionBoudary = CB(i,:);
        [A,Q] = BoundaryCondition(ConnectionBoudary,rho,beta,A0,A,Q,x,dt,L);
    end
    
    [A,Q] = MacCormack_CNMethod55(number_of_vessel,A,Q,beta,x,L,dt,A0,Rt,rho,dx,E,a,M,MA,K,t,N,init_A,init_Q,idxOutflowArtery);
    
    captureQ1(time_step) = Q{1}(20);
    captureQ2(time_step) = Q{8}(117);
    captureQ3(time_step) = Q{37}(5);
    captureQ4(time_step) = Q{54}(161);
    
    captureA1(time_step) = A{1}(20);
    captureA2(time_step) = A{8}(117);
    captureA3(time_step) = A{37}(5);
    captureA4(time_step) = A{54}(161);
end

%%
plot(x{1},Q{1})
