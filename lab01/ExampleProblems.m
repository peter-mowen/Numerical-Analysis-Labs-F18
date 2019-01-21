%% Example Problem 1
%plots the graph of cos(x) on the interval [0, 2*pi]

%creates a vector xx for plotting
xx = 0:0.01:2*pi;

%creates a function that can be evaluated at specific x values
f = @(x) cos(x);

%plots the graph of f(x) with a solid blue curve
plot(xx, f(xx), 'k','LineWidth', 2) 
%sets axis limits
axis([0, 2*pi, -1.5, 1.5])
legend({'cos(x)'}, 'location', 'SouthWest')

%% Example Problem 2
%chooses two random integers from [-10, 10] and checks the sign of the product

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

%some formats for fprintf: 
%  %d - whole number, 
%  %m.nf - floating-point number with m digits before decimal, n after, 
%  \n - carriage return

