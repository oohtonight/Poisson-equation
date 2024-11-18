close all; clc;

% Parameters
N = 1000;
h = 1 / N; 
x = linspace(0, 1, N+1)'; % Grid points
u0 = 0; % Boundary condition: u(0)
u1 = 0; % Boundary condition: u(1)

% Source function f(x)
f = x.^2; 

% Constructing the matrix A
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


b = f(2:end-1); 
b(1) = b(1) - u0 / h^2; 
b(end) = b(end) - u1 / h^2; 

% LU Decomposition
[L, U, P] = lu(A); % P is permutation matrix


y = L \ (P * b); 
u_interior = U \ y; 
u_numerical = [u0; u_interior; u1]; 

u_analytical = -x.^4 / 12 + x / 12;

error = norm(u_numerical - u_analytical, Inf);
fprintf('Error: %e\n', error);

figure;
plot(x, u_numerical, 'r-', 'LineWidth', 2); hold on;
plot(x, u_analytical, 'b--', 'LineWidth', 2);
xlabel('x', 'FontSize', 12);
ylabel('u(x)', 'FontSize', 12);
legend('Numerical Solution', 'Analytical Solution');
title('Case when f(x) = x^2', 'FontSize', 14);
grid on;
hold off;
