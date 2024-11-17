close all; clc;

% Parameters
N = 1000; % Number of intervals
h = 1 / N; % Grid step
x = linspace(0, 1, N+1)'; % Grid points
u0 = 0; % Boundary condition: u(0)
u1 = 0; % Boundary condition: u(1)

% Source function f(x)
f = x.^2; % f(x) = x^2

% Constructing the matrix
A = zeros(N-1, N-1);
for i = 1:N-1
    if i > 1
        A(i, i-1) = -1 / h^2; % Lower diagonal
    end
    A(i, i) = 2 / h^2; % Main diagonal
    if i < N-1
        A(i, i+1) = -1 / h^2; % Upper diagonal
    end
end

% Right-hand side b
b = f(2:end-1); % Right-hand side for internal points
b(1) = b(1) - u0 / h^2; % Adjusting for u(0)
b(end) = b(end) - u1 / h^2; % Adjusting for u(1)

% Calculating the numerical solution
u_interior = A \ b; % Solution for interior points
u_numerical = [u0; u_interior; u1]; % Complete numerical solution

% Updated analytical solution
u_analytical = -x.^4 / 12 + x / 12;

% Calculating the error
error = norm(u_numerical - u_analytical, Inf);
fprintf('Maximum error: %e\n', error);

% Plotting
figure;
plot(x, u_numerical, 'r-', 'LineWidth', 2); hold on;
plot(x, u_analytical, 'b--', 'LineWidth', 2);
xlabel('x', 'FontSize', 12);
ylabel('u(x)', 'FontSize', 12);
legend('Numerical Solution', 'Analytical Solution');
title('Case when f(x) = x^2', 'FontSize', 14);
grid on;
hold off;
