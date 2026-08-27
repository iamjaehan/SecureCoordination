function a = solveSchedule(tauHat, tauTilde, epsVec, lambda, sMin)
%SOLVESCHEDULE Robust min-max vertiport sequencing (paper Sec. IV-E, Eq. 12).
%
%   min_a max_{delta_i in [-eps_i,eps_i]}
%       sum_i (1-lambda_i)*(a_i - (tauTilde_i+delta_i))^2 + lambda_i*(a_i-tauHat_i)^2
%   s.t. a_i - a_{i-1} >= sMin,  i = 2,...,N
%
%   The inner maximization over delta_i has the closed form
%       max_{|delta_i|<=eps_i} (a_i-tauTilde_i-delta_i)^2 = (|a_i-tauTilde_i|+eps_i)^2,
%   so the problem reduces to a convex QP in (a,t) with epigraph variable
%       t_i >= |a_i - tauTilde_i|.
%
%   lambda(i) == 0 -> vehicle i fully distrusted (worst-case surveillance term only)
%   lambda(i) == 1 -> vehicle i fully trusted (reported ETA term only)

N = numel(tauHat);
tauHat   = tauHat(:);
tauTilde = tauTilde(:);
epsVec   = epsVec(:);
lambda   = lambda(:);

% variables x = [a(1:N); t(1:N)]
H = zeros(2*N);
f = zeros(2*N,1);
for i = 1:N
    H(i,i)     = 2*lambda(i);
    H(N+i,N+i) = 2*(1-lambda(i));
    f(i)       = -2*lambda(i)*tauHat(i);
    f(N+i)     = 2*(1-lambda(i))*epsVec(i);
end

rows = zeros(2*N + (N-1), 2*N);
rhs  = zeros(2*N + (N-1), 1);
r = 0;
for i = 1:N
    r = r+1; rows(r,i) = 1;  rows(r,N+i) = -1; rhs(r) =  tauTilde(i); % t_i >= a_i-tauTilde_i
    r = r+1; rows(r,i) = -1; rows(r,N+i) = -1; rhs(r) = -tauTilde(i); % t_i >= tauTilde_i-a_i
end
for i = 2:N
    r = r+1; rows(r,i-1) = 1; rows(r,i) = -1; rhs(r) = -sMin; % a_i-a_{i-1} >= sMin
end
A = rows;
b = rhs;

lb = [-inf(N,1); zeros(N,1)];
ub = inf(2*N,1);

opts = optimoptions('quadprog','Display','off');
x = quadprog(H, f, A, b, [], [], lb, ub, [], opts);
a = x(1:N);

end
