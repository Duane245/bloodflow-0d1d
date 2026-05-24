%%
% W1parentn1 = interp1(x{ConnectionBoudary(1)},W1.artery1,L(ConnectionBoudary(1))-lamda1.artery1(end)*dt);
% W2daughter1n1 = interp1(x{ConnectionBoudary(2)},W2.artery2,-lamda2.artery2(1)*dt);
% W2daughter2n1 = interp1(x{ConnectionBoudary(3)},W2.artery3,-lamda2.artery3(1)*dt);

str1 = ['W1parentn1 = interp1(x{ConnectionBoudary(1)},W1.artery' num2str(ConnectionBoudary(1)) ',L(ConnectionBoudary(1))-lamda1.artery' num2str(ConnectionBoudary(1)) '(end)*dt);'];
eval(str1);
str2 = ['W2daughter1n1 = interp1(x{ConnectionBoudary(2)},W2.artery' num2str(ConnectionBoudary(2)) ',-lamda2.artery' num2str(ConnectionBoudary(2)) '(1)*dt);'];
eval(str2);
str3 = ['W2daughter2n1 = interp1(x{ConnectionBoudary(3)},W2.artery' num2str(ConnectionBoudary(3)) ',-lamda2.artery' num2str(ConnectionBoudary(3)) '(1)*dt);'];
eval(str3);

F = @(X) [abs(X(2))-abs(X(4))-abs(X(6));
    0.5*rho*(X(2)/X(1))^2+abs(beta(ConnectionBoudary(1))*((sqrt(abs(X(1))))-(sqrt(A0(ConnectionBoudary(1)))))*1e4)-0.5*rho*(X(4)/X(3))^2-abs(beta(ConnectionBoudary(2))*((sqrt(abs(X(3))))-(sqrt(A0(ConnectionBoudary(2)))))*1e4);
    0.5*rho*(X(2)/X(1))^2+abs(beta(ConnectionBoudary(1))*((sqrt(abs(X(1))))-(sqrt(A0(ConnectionBoudary(1)))))*1e4)-0.5*rho*(X(6)/X(5))^2-abs(beta(ConnectionBoudary(3))*((sqrt(abs(X(5))))-(sqrt(A0(ConnectionBoudary(3)))))*1e4);
    -abs(X(2)/X(1))-4*sqrt(beta(ConnectionBoudary(1))/2/rho*sqrt(abs(X(1))))*100+W1parentn1;
    -abs(X(4)/X(3))+4*sqrt(beta(ConnectionBoudary(2))/2/rho*sqrt(abs(X(3))))*100+W2daughter1n1;
    -abs(X(6)/X(5))+4*sqrt(beta(ConnectionBoudary(3))/2/rho*sqrt(abs(X(5))))*100+W2daughter2n1];
% X0 = [A.artery1(end);Q.artery1(end);A.artery2(1);Q.artery2(1);A.artery3(1);Q.artery3(1)];

str = ['X0 = [A.artery' num2str(ConnectionBoudary(1)) '(end);Q.artery' num2str(ConnectionBoudary(1)) '(end);A.artery' num2str(ConnectionBoudary(2)) '(1);Q.artery' num2str(ConnectionBoudary(2)) '(1);A.artery' num2str(ConnectionBoudary(3)) '(1);Q.artery' num2str(ConnectionBoudary(3)) '(1)];'];
eval(str);

options = optimoptions('fsolve','Display','iter');
% options = optimoptions('fsolve','Display','none');
[X,fval] = fsolve(F,X0,options);

% A.artery1(end) = abs(X(1));
% Q.artery1(end) = abs(X(2));
% A.artery2(1) = abs(X(3));
% Q.artery2(1) = abs(X(4));
% A.artery3(1) = abs(X(5));
% Q.artery3(1) = abs(X(6));

str1 = ['A.artery' num2str(ConnectionBoudary(1)) '(end) = abs(X(1));'];
eval(str1);
str2 = ['Q.artery' num2str(ConnectionBoudary(1)) '(end) = abs(X(2));'];
eval(str2);
str3 = ['A.artery' num2str(ConnectionBoudary(2)) '(1) = abs(X(3));'];
eval(str3);
str4 = ['Q.artery' num2str(ConnectionBoudary(2)) '(1) = abs(X(4));'];
eval(str4);
str5 = ['A.artery' num2str(ConnectionBoudary(3)) '(1) = abs(X(5));'];
eval(str5);
str6 = ['Q.artery' num2str(ConnectionBoudary(3)) '(1) = abs(X(6));'];
eval(str6);

% Parent outflow boundary condition
% W1.artery1(N(ConnectionBoudary(1))) = Q.artery1(end)./A.artery1(end) + 4*sqrt(beta(ConnectionBoudary(1))/2/rho*sqrt(abs(A.artery1(end))))*100;
% W2.artery1(N(ConnectionBoudary(1))) = Q.artery1(end)./A.artery1(end) - 4*sqrt(beta(ConnectionBoudary(1))/2/rho*sqrt(abs(A.artery1(end))))*100;
str1 = ['W1.artery' num2str(ConnectionBoudary(1)) '(N(ConnectionBoudary(1))) = Q.artery' num2str(ConnectionBoudary(1)) '(end)./A.artery' num2str(ConnectionBoudary(1)) '(end) + 4*sqrt(beta(ConnectionBoudary(1))/2/rho*sqrt(abs(A.artery' num2str(ConnectionBoudary(1)) '(end))))*100;'];
str2 = ['W2.artery' num2str(ConnectionBoudary(1)) '(N(ConnectionBoudary(1))) = Q.artery' num2str(ConnectionBoudary(1)) '(end)./A.artery' num2str(ConnectionBoudary(1)) '(end) - 4*sqrt(beta(ConnectionBoudary(1))/2/rho*sqrt(abs(A.artery' num2str(ConnectionBoudary(1)) '(end))))*100;'];
eval(str1);
eval(str2);

% daughter1 inflow boundary condition (n time step W2 computing n+1 timestep W1)
% W1.artery2(1) = Q.artery2(1)./A.artery2(1) + 4*sqrt(beta(ConnectionBoudary(2))/2/rho*sqrt(abs(A.artery2(1))))*100;
% W2.artery2(1) = Q.artery2(1)./A.artery2(1) - 4*sqrt(beta(ConnectionBoudary(2))/2/rho*sqrt(abs(A.artery2(1))))*100;
str1 = ['W1.artery' num2str(ConnectionBoudary(2)) '(1) = Q.artery' num2str(ConnectionBoudary(2)) '(1)./A.artery' num2str(ConnectionBoudary(2)) '(1) + 4*sqrt(beta(ConnectionBoudary(2))/2/rho*sqrt(abs(A.artery' num2str(ConnectionBoudary(2)) '(1))))*100;'];
str2 = ['W2.artery' num2str(ConnectionBoudary(2)) '(1) = Q.artery' num2str(ConnectionBoudary(2)) '(1)./A.artery' num2str(ConnectionBoudary(2)) '(1) - 4*sqrt(beta(ConnectionBoudary(2))/2/rho*sqrt(abs(A.artery' num2str(ConnectionBoudary(2)) '(1))))*100;'];
eval(str1);
eval(str2);

% daughter2 inflow boundary condition (n time step W2 computing n+1 timestep W1)
% W1.artery3(1) = Q.artery3(1)./A.artery3(1) + 4*sqrt(beta(ConnectionBoudary(3))/2/rho*sqrt(abs(A.artery3(1))))*100;
% W2.artery3(1) = Q.artery3(1)./A.artery3(1) - 4*sqrt(beta(ConnectionBoudary(3))/2/rho*sqrt(abs(A.artery3(1))))*100;

str1 = ['W1.artery' num2str(ConnectionBoudary(3)) '(1) = Q.artery' num2str(ConnectionBoudary(3)) '(1)./A.artery' num2str(ConnectionBoudary(3)) '(1) + 4*sqrt(beta(ConnectionBoudary(3))/2/rho*sqrt(abs(A.artery' num2str(ConnectionBoudary(3)) '(1))))*100;'];
str2 = ['W2.artery' num2str(ConnectionBoudary(3)) '(1) = Q.artery' num2str(ConnectionBoudary(3)) '(1)./A.artery' num2str(ConnectionBoudary(3)) '(1) - 4*sqrt(beta(ConnectionBoudary(3))/2/rho*sqrt(abs(A.artery' num2str(ConnectionBoudary(3)) '(1))))*100;'];
eval(str1);
eval(str2);
