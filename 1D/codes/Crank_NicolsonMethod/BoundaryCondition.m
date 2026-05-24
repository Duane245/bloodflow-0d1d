function [A,Q] = BoundaryCondition(ConnectionBoudary,rho,beta,A0,A,Q,x,dt,L)
%BOUNDARYCONDITION  3-vessel junction (parent + 2 daughters) coupling.
%   Solves 6 nonlinear equations:
%     mass conservation (Qp = Qd1 + Qd2)
%     total-pressure continuity across the junction (2 eqs)
%     outgoing characteristic compatibility (3 eqs)
%   Newton with analytical Jacobian + inline 2-pt linear interp.

ip  = ConnectionBoudary(1);
id1 = ConnectionBoudary(2);
id2 = ConnectionBoudary(3);

bp  = beta(ip);  bd1 = beta(id1); bd2 = beta(id2);
A0p = A0(ip);    A0d1= A0(id1);   A0d2= A0(id2);
sA0p  = sqrt(A0p);
sA0d1 = sqrt(A0d1);
sA0d2 = sqrt(A0d2);

% characteristic speeds at junction
lamda1end = Q{ip}(end)/A{ip}(end)  + sqrt(bp /(2*rho)*sqrt(A{ip}(end))) *100;
lamda20   = Q{id1}(1) /A{id1}(1)   - sqrt(bd1/(2*rho)*sqrt(A{id1}(1))) *100;
lamda30   = Q{id2}(1) /A{id2}(1)   - sqrt(bd2/(2*rho)*sqrt(A{id2}(1))) *100;

% foot-of-characteristic interpolations (inline 2-pt linear on uniform grid)
xp  = x{ip};   Ap_arr  = A{ip};   Qp_arr  = Q{ip};
xd1 = x{id1};  Ad1_arr = A{id1};  Qd1_arr = Q{id1};
xd2 = x{id2};  Ad2_arr = A{id2};  Qd2_arr = Q{id2};

dxp  = xp(2);  np  = numel(xp);
dxd1 = xd1(2); nd1 = numel(xd1);
dxd2 = xd2(2); nd2 = numel(xd2);

xq = L(ip) - lamda1end*dt;
pos = xq/dxp + 1;  i0 = floor(pos); if i0<1, i0=1; elseif i0>np-1, i0=np-1; end
fr = pos - i0;
Ap_foot = (1-fr)*Ap_arr(i0) + fr*Ap_arr(i0+1);
Qp_foot = (1-fr)*Qp_arr(i0) + fr*Qp_arr(i0+1);

xq = -lamda20*dt;
pos = xq/dxd1 + 1; i0 = floor(pos); if i0<1, i0=1; elseif i0>nd1-1, i0=nd1-1; end
fr = pos - i0;
Ad1_foot = (1-fr)*Ad1_arr(i0) + fr*Ad1_arr(i0+1);
Qd1_foot = (1-fr)*Qd1_arr(i0) + fr*Qd1_arr(i0+1);

xq = -lamda30*dt;
pos = xq/dxd2 + 1; i0 = floor(pos); if i0<1, i0=1; elseif i0>nd2-1, i0=nd2-1; end
fr = pos - i0;
Ad2_foot = (1-fr)*Ad2_arr(i0) + fr*Ad2_arr(i0+1);
Qd2_foot = (1-fr)*Qd2_arr(i0) + fr*Qd2_arr(i0+1);

W1p  = Qp_foot /Ap_foot  + 4*sqrt(bp /(2*rho)*sqrt(Ap_foot ))*100;
W2d1 = Qd1_foot/Ad1_foot - 4*sqrt(bd1/(2*rho)*sqrt(Ad1_foot))*100;
W2d2 = Qd2_foot/Ad2_foot - 4*sqrt(bd2/(2*rho)*sqrt(Ad2_foot))*100;

% Constants for the characteristic terms 4*k*X^(1/4)
kp  = sqrt(bp /(2*rho))*100;
kd1 = sqrt(bd1/(2*rho))*100;
kd2 = sqrt(bd2/(2*rho))*100;

% Newton on X = [Ap, Qp, Ad1, Qd1, Ad2, Qd2]
X = [A{ip}(end); Q{ip}(end); A{id1}(1); Q{id1}(1); A{id2}(1); Q{id2}(1)];

J = zeros(6,6);
J(1,2) = 1; J(1,4) = -1; J(1,6) = -1;

tol = 1e-10;
for iter = 1:20
    X1 = X(1); X2 = X(2); X3 = X(3); X4 = X(4); X5 = X(5); X6 = X(6);
    sX1 = sqrt(X1); sX3 = sqrt(X3); sX5 = sqrt(X5);

    F = [X2 - X4 - X6;
         0.5*rho*(X2/X1)^2 + bp *(sX1 - sA0p )*1e4 - 0.5*rho*(X4/X3)^2 - bd1*(sX3 - sA0d1)*1e4;
         0.5*rho*(X2/X1)^2 + bp *(sX1 - sA0p )*1e4 - 0.5*rho*(X6/X5)^2 - bd2*(sX5 - sA0d2)*1e4;
         -X2/X1 - 4*kp *X1^0.25 + W1p;
         -X4/X3 + 4*kd1*X3^0.25 + W2d1;
         -X6/X5 + 4*kd2*X5^0.25 + W2d2];

    if max(abs(F)) < tol, break; end

    J(2,1) = -rho*X2^2/X1^3 + 0.5*bp /sX1*1e4;
    J(2,2) =  rho*X2/X1^2;
    J(2,3) =  rho*X4^2/X3^3 - 0.5*bd1/sX3*1e4;
    J(2,4) = -rho*X4/X3^2;

    J(3,1) = J(2,1);
    J(3,2) = J(2,2);
    J(3,5) =  rho*X6^2/X5^3 - 0.5*bd2/sX5*1e4;
    J(3,6) = -rho*X6/X5^2;

    J(4,1) =  X2/X1^2 - kp *X1^(-0.75);
    J(4,2) = -1/X1;

    J(5,3) =  X4/X3^2 + kd1*X3^(-0.75);
    J(5,4) = -1/X3;

    J(6,5) =  X6/X5^2 + kd2*X5^(-0.75);
    J(6,6) = -1/X5;

    X = X - J\F;
end

A{ip}(end) = X(1); Q{ip}(end) = X(2);
A{id1}(1)  = X(3); Q{id1}(1)  = X(4);
A{id2}(1)  = X(5); Q{id2}(1)  = X(6);
end
