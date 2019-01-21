% NOTE: The bisection method is implemented in the BisectionMethod.m
%       function file. Please look there for implementation of Bisection method!
%% Problem 1 
%{
Use the Bisection method to find a solution accurate to within epsilon = 10^(-8) 
    for x - 2^(-x) = 0 on the interval [0, 1]. 
    Set the maximum number of iterations to be 30, and use the
    stopping criteria (b_n - a_n)/2 < epsilon
%}

fprintf('Problem 1\n')
% Set parameters for Bisection Method input
lowerBound = 0;         % set lower bound
upperBound = 1;         % set upper bound
tolerance = 10^(-8);    % set tolerance
stopCriteria = 1;       % see case 1 of stopCriteria in BisectionMethod.m
maxIterations = 30;     % set maximum iterations
f = @(x) x - 2^(-x);    % define function f

% Run bisection method
[solved, approxSolution, numIterations] = ...
    BisectionMethod( lowerBound, upperBound, tolerance, stopCriteria, maxIterations, f);

% Print Results
if solved == 0
    fprintf('Method could not find a zero within %d iterations for:\n', maxIterations)
    fprintf('\tepsilon = %s\n', tolerance)
else
    fprintf('\tapproximation: %1.8f \n\titerations: %d \n\ttolerance: %s\n',approxSolution, numIterations, tolerance)
end
%% Problem 2
%{
Repeat Problem 1 using epsilon = 10^-12.
%}

fprintf('\nProblem 2\n')

% Set parameters for Bisection Method input
lowerBound = 0;         % set lower bound
upperBound = 1;         % set upper bound
tolerance = 10^(-12);   % set tolerance
stopCriteria = 1;       % see case 1 of stopCriteria in BisectionMethod.m
maxIterations = 30;     % set maximum iterations
f = @(x) x - 2^(-x);    % define function f

% Run Bisection Method
[solved, approxSolution, numIterations] = ...
    BisectionMethod(lowerBound, upperBound, tolerance, stopCriteria, maxIterations, f);

% Print Results
if solved == 0
    fprintf('Method could not find a zero within %d iterations for:\n', maxIterations)
    fprintf('\tepsilon = %s\n', tolerance)
else
    fprintf('\tapproximation: %1.8f \n\titerations: %d \n\ttolerance: %f\n',approxSolution, numIterations, tolerance)
end
%% Problem 3
%{   
Plot the graphs of y = x and y = 2*sin(x). Use the Bisection method to find an approximation
    to within epsilon = 10^-8 to the first positive value of x with x = 2sin x. 
    Use the stopping criteria: |(pn - pnMinus1)/pn| < epsilon

Remark: The first positive value of x = 2sin(x) will occur on [0,2[ because
        2sin(x) < 2 for every x
%}

fprintf('\nProblem 3\n')

% Define functions as 
iota = @(x) x;                  % define function iota as the identity function (y=x)
alpha = @(x) 2*sin(x);          % define function alpha

% Define set of inputs and outpoints for functions
xValues = 0:0.001:2*pi;         % define inputs for both functions iota and alpha
yValuesIota = xValues;          % because y = x for iota
yValuesAlpha = alpha(xValues);  % find y values for alpha

% Plot functions
plot(xValues,yValuesAlpha)      % plot alpa
axis([0, 2*pi, -2.5, 2.5])      % set axis limits
hold on                         % MATLAB magic to print another function on same plots
plot(xValues, yValuesIota)      % plot iota on same graph as alpha
legend('alpha(x) = 2sinx', 'iota(x) = x'); xlabel('x'); ylabel('y')

% Set parameters for Bisection Method input
lowerBound = 1.5;               % set lower bound based on looking at graph
upperBound = 2;                 % set upper bound based on looking at graph
tolerance = 10^(-8);            % define tolerance
stopCriteria = 2;               % see case 2 of stopCriteria in BisectionMethod.m
maxIterations = 30;             % set max iterations
f = @(x) iota(x) - alpha(x);    % define function with fixed point at some p

%Run Bisection Method
[solved, approxSolution, numIterations] = ...
    BisectionMethod(lowerBound, upperBound, tolerance, stopCriteria, maxIterations, f);

% Print Results
if solved == 0
    fprintf('Method could not find a zero within %d iterations for:\n', maxIterations)
    fprintf('\tepsilon = %f\n', tolerance)
else
    fprintf('\tapproximation: %1.8f \n\titerations: %d \n\ttolerance: %s\n',approxSolution, numIterations, tolerance)
end
%% Probelm 4
%{
Find an approximation to the 3root(25) correct to within epsilon = 10^(-10) 
    using the Bisection method. Hint: Consider f(x) = x^3 - 25. Use the 
    stopping criteria:
        abs((pn - pnMinus1)/pn) < epsilon

Remark: we can define a function with a fixed point at root(3,25) by
        noting that if
            f := x -> root(3,x) : R -> R and   f(25) = root(3,25)
        => f^3:= x ->      x    : R -> R and f^3(25) = 25 
           and if g := x -> f^3(x) - 25: R -> R  and if g(x) = 0
        => f^3(x) - 25 = 0, which is a root finding problem whose solution
           is the answer for which we are searching and can be solved with
           the bisection method using this definition of g.

Remark:    8 < 25 < 27 
        => root(3,8) < root(3,25) < root(3,27) 
        => 2 < root(3,25) < 3
        => root(3,25) in ]2,3[
%}

fprintf('\nProblem 4\n')

lowerBound = 2;                 % set lower bound based on looking at graph
upperBound = 3;                 % set upper bound based on looking at graph
tolerance = 10^(-8);            % set tolerance
stopCriteria = 2;               % see case 2 of stopCriteria in BisectionMethod.m
maxIterations = 30;             % set max iterations
g = @(x) x^3 - 25;              % define function with fixed point at some p

[solved, approxSolution, numIterations] = ...
    BisectionMethod(lowerBound, upperBound, tolerance, stopCriteria, maxIterations, g);
% Print Results
if solved == 0
    fprintf('Method could not find a zero within %d iterations for:\n', maxIterations)
    fprintf('\tepsilon = %s\n', tolerance)
else
    fprintf('\tapproximation: %1.8f \n\titerations: %d \n\ttolerance: %s\n',approxSolution, numIterations, tolerance)
end