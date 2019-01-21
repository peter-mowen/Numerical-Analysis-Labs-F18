%% MATLAB in a nutshell
% Note this code is to be stored in a file: nutshell.m
% This is a very brief and quick introduction to MATLAB
% Read also the corresponding sections in the MATLAB documentation
% MATLAB documentation can be found here: https://www.mathworks.com/help/matlab

%%% variables
alpha1 = 10 % scalar
vec1 = [1, 2, 3] % row vector
vec2 = [1; 2; 3] % column vector
xValues = 0:0.01:2*pi % range vector

%%% semicolon and comments
a = 10; % will define a and not print the value of a
a = 10 % will defien a and print the value of a
% everything after the percentage sign will be ignored
%{
 This is a multiline comment
 available from MATLAB 7 (R14)
 and may help to comment out a
 larger block of code instead of
 using single percentage signs
%}

%%% functions and evaluation
f = @(x) cos(x); % defines f(x) = cos(x)
f(2) % evaluates f(2), i.e. cos(2) 
yValues = f(xValues); % evaluates all xValues

% 2d plotting
plot(xValues, yValues, 'k', 'LineWidth', 2) % What is 'k'?
axis([0, 2*pi, -1.5, 1.5]) % set axis limits
legend({'cos(x)'}, 'location', 'SouthWest') % legend position

%%% getting help / documentation
help plot
help doc 

%%% conditional statements
a = randi([-10 10]);
b = randi([-10 10]);

fprintf('a = %d, b = %d, a*b = %d\n', a, b, a*b)

if a*b == 0
    fprintf('The product of a and b is zero\n')
elseif a*b <0
    fprintf('The product of a and b is negative\n')
else
    fprintf('The product of a and b is positive\n')
end

%%% for loop
for n = 1:1:10
   factorial(n)
end

%%% while loop
n = 1
while n <= 10
  factorial(n)
  n = n + 1;
end

%%% fprintf statement (\n: carriage return)
fprintf('This is some text\n'); % print the plain text
fprintf('A random number: %d\n', randi([0, 10])); % print a random whole number between 0 and 10
format long % show 16 significant digits 
fprintf('Pi is: %f\n', pi);
fprintf('Pi is: %1.2f\n', pi) % 1 digit before decimal and 2 after
format short % show 6 significant digits
fprintf('Pi: %f\n', pi);

%%% publish this MATLAB file (default HTML)
publish('class1.m','pdf') % PDF
publish('class1.m','doc') % MS Word
