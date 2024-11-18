close all; clc;

% Parameters
N = 1000; 
h = 1 / N; 
x = linspace(0, 1, N+1)'; 
u0 = 0; % Boundary condition: u(0)
u1 = 0; % Boundary condition: u(1)

% Source function f(x)
f = sin(pi * x); 

% Constructing a fourth-order matrix
A = zeros(N-1, N-1);
for i = 1:N-1
    if i > 2
        A(i, i-2) = -1 / (12 * h^2); 
    end
    if i > 1
        A(i, i-1) = 16 / (12 * h^2); 
    end
    A(i, i) = -30 / (12 * h^2); 
    if i < N-1
        A(i, i+1) = 16 / (12 * h^2); 
    end
    if i < N-2
        A(i, i+2) = -1 / (12 * h^2); 
    end
end


b = f(2:end-1); 
b(1) = b(1) - u0 * (16 / (12 * h^2)); 
b(end) = b(end) - u1 * (16 / (12 * h^2)); 

% LU Decomposition
[L, U, P] = lu(A); 


y = L \ (P * b); 
u_interior = U \ y; 
u_numerical = [u0; u_interior; u1]; 


u_analytical = -sin(pi * x) / pi^2;

error = norm(u_numerical - u_analytical, Inf);
fprintf('Error: %e\n', error);


figure;
plot(x, u_numerical, 'r-', 'LineWidth', 2); hold on;
plot(x, u_analytical, 'b--', 'LineWidth', 2);
xlabel('x', 'FontSize', 12);
ylabel('u(x)', 'FontSize', 12);
legend('Numerical Solution (Fourth-Order)', 'Analytical Solution');
title('Case when f(x) = sin(\pi x)', 'FontSize', 14);
grid on;
hold off;
