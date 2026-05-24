function [W1,W2] = CorStep(W1i,W2i,W1,W2,lamda1,lamda2,Q,A,Cf,Cv,dx,dt,N,A0)
dQ2(2:N-1) = (Q(3:N)-2*Q(2:N-1)+Q(1:N-2))/(dx^2);

W1(2:N-1) = (W1(2:N-1)+W1i(2:N-1))/2 - (dt/dx/2)*(W1(2:N-1)-W1(1:N-2)).*lamda1(2:N-1) - (dt/2)*(Cf/A0)*Q(2:N-1)./A(2:N-1) + dt/2*Cv*dQ2(2:N-1)./A(2:N-1);
W2(2:N-1) = (W2(2:N-1)+W2i(2:N-1))/2 - (dt/dx/2)*(W2(2:N-1)-W2(1:N-2)).*lamda2(2:N-1) - (dt/2)*(Cf/A0)*Q(2:N-1)./A(2:N-1) + dt/2*Cv*dQ2(2:N-1)./A(2:N-1);

end
