function [found, approxSolution, numIterations, error] = ...
    BisectionMethod(lowerBound,upperBound, tolerance, stopCriteriaCase, ...
    maxIterations, f)
%This function implements finding zeros using the Bisection Method
%   Inputs: lowerBound, upperBound, tolerance, stop criteria case number (1,2),
%           maxIterations, and a function f.
%   Output: found: a boolean value telling whether or not a root was found, 
%           approximate solution, number of iterations to reach solution
%   DEFINITION: Stoping criteria functions defined as follows: 
%       1) (upperBound - lowerBound)/2
%       2) |(pn - pnMinus1)/pn|

approxSolution = '';    % intialized incase of error
numIterations = '';     % intialized incase of error
found = 0;              % intialize found to false

% check to see if stopCriteriaCase is valid. Return if it is not.
if stopCriteriaCase == 1
    stopCriteriaF = @(lowerBound,upperBound) (upperBound - lowerBound)/2;
elseif stopCriteriaCase == 2
    stopCriteriaF = @(pn, pnMinus1) abs((pn - pnMinus1)/pn);
else
    fprintf('Invalid stop criteria case.\n')
    return
end
n = 1; % initialize counter
error = 10^8; % junk value to intialize error variable
while n <= maxIterations
    % find midpoint
    pn = (lowerBound + upperBound)/2;
    %test f(p) and stopping criteria.
    if stopCriteriaCase == 1
        error = stopCriteriaF(lowerBound,upperBound);
        if f(pn) == 0 || error < tolerance
            found = 1;
            approxSolution = pn;
            numIterations = n;
            return
        end
    elseif stopCriteriaCase == 2 && n > 1
        error = stopCriteriaF(pn,pnMinus1);
        if f(pn) == 0 || error < tolerance
            found = 1;
            approxSolution = pn;
            numIterations = n;
            return
        end
    end
    
    n = n + 1; % increment n
    fprintf('% d, %.8f, %.8f, %.8f, %.8f\n', n, lowerBound, upperBound, pn, error);
    % assign new bounderies and save previous p
    sgn = sign(f(lowerBound))*sign(f(pn));
    if sgn < 0
        pnMinus1 = upperBound; % used in stopping criteria 2
        upperBound = pn; % update upperBound
    else
        pnMinus1 = lowerBound; % used in stopping criteria 2
        lowerBound = pn; % update lowerBound
    end
end