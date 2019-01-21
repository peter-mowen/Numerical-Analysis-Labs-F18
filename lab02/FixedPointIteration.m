function [foundSol, numIterations, approxSol, finalError] = FixedPointIteration(p0, tolerance,...
    stoppingCriteria, maxIterations, f)
% This function finds an approximate zero to the input function using
%   fixed-point iteration.
%   INPUTS: po - starting guess, tolerance, maximum number of iterations,
%           a function f
%   OUTPUTS: foundSol - a boolean value of wether or not a solution was 
%           found within given restrainst, approxSol - approximate solution 
%           of input function near p0.
foundSol = 0; % we haven't found a solution yet
n = 0; % initialize iteration counter
pn = p0; %set pn to p0 for first run
while n < maxIterations
    xn = f(pn);
    if stoppingCriteria == 1
        error = absoluteError(xn, pn);
        if error < tolerance
            foundSol = 1; % solution found
            numIterations = n;
            approxSol = xn; % appoximate solution
            finalError = error;
            break
        end
    end
    n = n + 1; % increment iteration counter
    pn = xn; % set input for next iteration to output from this one
    fprintf('n: %d \tp%d: %.10f \t |error|: %.10f\n', n, n, xn, error);
end
end