%% Problem 1
% For the linear system:
clear
clc
A = [ 10, -1,  2,  0;
      -1, 11, -1,  3;
       2, -1, 10, -1;
       0,  3, -1.  8];

b = [6, 25, -11, 15]';
x0 = zeros(4,1);
tolerance = 10^-3;

%(a) use the Jacobi method to solve the linear system
[solJacobi, nJacobi] = JacobisMethod(A, b, x0, tolerance)

%(b) repeat part (a) using the Gauss-Seidel method
[solGS, nGS] = GSMethod(A, b, x0, tolerance)

solExact = A\b

%% Problem 2
% For the linear system


A = [ 4,  -1,  0,  0,   0,  0;
      -1,  4, -1,  0,   0,  0;
       0, -1,  4,  0,   0,  0;
       0,  0,  0,  4,  -1,  0;
       0,  0,  0, -1,   4, -1;
       0,  0,  0,  0,  -1,  4];

b = [0, 5, 0, 6, -2, 6]';
x0 = zeros(6,1);
tolerance = 10^-4;

%(a) use the Jacobi method to solve the linear system
[solJacobi, nJacobi] = JacobisMethod(A, b, x0, tolerance)

%(b) repeat part (a) using the Gauss-Seidel method
[solGS, nGS] = GSMethod(A, b, x0, tolerance)

solExact = A\b

%% Problem 3

n = 80;
A = zeros(n);

for i = 1: n
    for j = 1: n
        if i == j
            A(i,j) = i;
        elseif ( j == i + 2 ) && ( 1 <= i ) && ( i <= n-2 )
            A(i,j) = 0.5 * i;
        elseif ( j == i - 2 ) && ( 3 <= i ) && ( i <= n   )
            A(i,j) = 0.5 * i;
        elseif ( j == i + 4 ) && ( 1 <= i ) && ( i <= n-4 )
            A(i,j) == 0.25*i;
        elseif ( j == i - 4 ) && ( 5 <= i ) && ( i <= n   )
            A(i,j) == 0.25*i;
        else
            A(i,j) = 0;
        end
    end
end

x0 = zeros(n,1);

b = x0 + pi;

tolerance = 1e-5;

[sol, n, Tj] = JacobisMethod(A, b, x0, tolerance);
sol

[sol, n, Tgs] = GSMethod(A, b, x0, tolerance);
sol

%% Problem 4
clc
A = [ 2, -1, 1;
      2,  2, 2;
     -1, -1, 2];

b = [1, 2, -1]';

x0 = zeros(3,1);

tolerance = 1e-5;

[sol, n, Tj] = JacobisMethod(A, b, x0, tolerance);
egnVlsTj = eig(Tj);
rhoOfTj = max(abs(egnVlsTj))

[sol, n, Tgs] = GSMethod(A, b, x0, tolerance);
egnVlsTgs = eig(Tgs);
rhoOfTj = max(abs(egnVlsTgs))

%% Problem 5
clc
A = [ 1, 2, -2;
      1, 1,  1;
      2, 2,  1];

b = [7, 2, 5]';

x0 = zeros(3,1);

tolerance = 1e-5;

[sol, n, Tj] = JacobisMethod(A, b, x0, tolerance);
egnVlsTj = eig(Tj);
rhoOfTj = max(abs(egnVlsTj))

[sol, n, Tgs] = GSMethod(A, b, x0, tolerance);
egnVlsTgs = eig(Tgs);
rhoOfTgs = max(abs(egnVlsTgs))

%% Jacobi's Method
function [sol, n, Tj] = JacobisMethod(A, b, x0, tolerance)
N = diag(diag(A));
P = N - A;

xn = x0;
error = 999999;
n = 1; % number of iterations
Tj = inv(N)*P;
c = inv(N)*b;
while error > tolerance; % TODO update to for loop with maxIter
xnPlus1 = Tj*xn + c;
error = CheckTolerance(xn, xnPlus1);
n = n+1;
xn = xnPlus1;
sol = xnPlus1;
end
end

%% Gauss-Seidel Method
function [sol, n, Tgs] = GSMethod(A, b, x0, tolerance)
D = diag(diag(A));
L = -(tril(A) - D);
U = -(triu(A) - D);

Tgs = inv(D-L)*U;
c = inv(D-L)*b;

xn = x0;
error = 999999;
n = 1; % number of iterations

while error > tolerance % TODO update to for loop with maxIter
xnPlus1 = Tgs*xn + c;
error = CheckTolerance(xn, xnPlus1);
n = n+1;
xn = xnPlus1;
end

sol = xnPlus1;

end

%% Check Tolerance Function
function error = CheckTolerance(xn, xnPlus1)
error = norm( xnPlus1 - xn, Inf)/norm(xnPlus1, Inf);
end