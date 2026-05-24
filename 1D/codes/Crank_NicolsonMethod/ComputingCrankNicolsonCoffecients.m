function [E,a,M,K,MA] = ComputingCrankNicolsonCoffecients(number_of_vessel,dt,dx,Cv,c0,N)
%COMPUTINGCRANKNICOLSONCOFFECIENTS  Pre-compute the per-vessel CN matrices.
%   M{i} is a tridiagonal (N-2)x(N-2) system matrix; MA{i} is the (N-2)xN
%   first-difference-over-two-cells operator. Both stored as sparse so
%   M{i}\K{i} uses the banded solver (Thomas) instead of dense LU.

E  = zeros(number_of_vessel,1);
a  = zeros(number_of_vessel,1);
K  = cell(number_of_vessel,1);
M  = cell(number_of_vessel,1);
MA = cell(number_of_vessel,1);

for i = 1:number_of_vessel
    E(i) = dt*Cv(i)/(dx^2);
    a(i) = -E(i)/2;
    b    = 1+E(i);
    c    = -dt*c0(i)^2/4/dx;

    n = N(i) - 2;

    % Tridiagonal M: diag = b, off-diag = a
    main = b * ones(n,1);
    off  = a(i) * ones(n,1);
    M{i} = spdiags([off, main, off], -1:1, n, n);

    K{i} = zeros(n,1);

    % MA is (n) x (N), with MA(k,k) = -c and MA(k,k+2) = c for k=1..n
    Nv   = N(i);
    rows = repmat((1:n).',1,2);
    cols = [(1:n).', (3:Nv).'];
    vals = repmat([-c, c], n, 1);
    MA{i} = sparse(rows(:), cols(:), vals(:), n, Nv);
end
end
