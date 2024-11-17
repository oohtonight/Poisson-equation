close all; clc;

% Parameters
N = 1000; % Number of intervals
h = 1 / N; % Step size
x = linspace(0, 1, N+1)'; % Grid points
u0 = 0; % Boundary condition: u(0)
u1 = 0; % Boundary condition: u(1)

% Source function f(x)
f = sin(pi * x); % f(x) = sin(pi * x)

% Constructing a fourth-order matrix
A = zeros(N-1, N-1);
for i = 1:N-1
    if i > 2
        A(i, i-2) = -1 / (12 * h^2); % -u_{i-2}
    end
    if i > 1
        A(i, i-1) = 16 / (12 * h^2); % 16u_{i-1}
    end
    A(i, i) = -30 / (12 * h^2); % -30u_i
    if i < N-1
        A(i, i+1) = 16 / (12 * h^2); % 16u_{i+1}
    end
    if i < N-2
        A(i, i+2) = -1 / (12 * h^2); % -u_{i+2}
    end
end

% Right-hand side b
b = f(2:end-1); % Right-hand side for interior points
b(1) = b(1) - u0 * (16 / (12 * h^2)); % Adjusting for u(0)
b(end) = b(end) - u1 * (16 / (12 * h^2)); % Adjusting for u(1)

% Calculating the numerical solution
u_interior = A \ b; % Solution for interior points
u_numerical = [u0; u_interior; u1]; % Complete numerical solution

% Analytical solution
u_analytical = -sin(pi * x) / pi^2;

% Calculating the error
error = norm(u_numerical - u_analytical, Inf);
fprintf('Maximum error: %e\n', error);

% Plotting
figure;
plot(x, u_numerical, 'r-', 'LineWidth', 2); hold on;
plot(x, u_analytical, 'b--', 'LineWidth', 2);
xlabel('x', 'FontSize', 12);
ylabel('u(x)', 'FontSize', 12);
legend('Numerical Solution (Fourth-Order)', 'Analytical Solution');
title('Case when f(x) = sin(\pi x)', 'FontSize', 14);
grid on;
hold off;
