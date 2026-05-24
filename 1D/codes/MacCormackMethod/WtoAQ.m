function [A,Q] = WtoAQ(W1,W2,rho,beta)
A = (W1-W2).^4/1024*(rho/beta)^2*1e-8;
Q = A.*(W1+W2)/2;
end
