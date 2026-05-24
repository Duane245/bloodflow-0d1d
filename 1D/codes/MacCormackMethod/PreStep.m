function [W1,W2] = PreStep(W1,W2,lamda1,lamda2,Q,A,Cf,Cv,dx,dt,N,A0)
dQ2(2:N-1) = (Q(3:N)-2*Q(2:N-1)+Q(1:N-2))/(dx^2);

W1(2:N-1) = W1(2:N-1) - (dt/dx)*(W1(3:N)-W1(2:N-1)).*lamda1(2:N-1) - dt*(Cf/A0)*Q(2:N-1)./A(2:N-1) + (dt*Cv)*dQ2(2:N-1)./A(2:N-1);
W2(2:N-1) = W2(2:N-1) - (dt/dx)*(W2(3:N)-W2(2:N-1)).*lamda2(2:N-1) - dt*(Cf/A0)*Q(2:N-1)./A(2:N-1) + (dt*Cv)*dQ2(2:N-1)./A(2:N-1);

end
