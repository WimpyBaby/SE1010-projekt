
% Initialize global variables

global c A v p a1 a2 g m L R df db h h1 b1 dh rd rb bb bd D d a gamma rh

c = 0.3;
A = 0.45;
v = 120/3.6;
p = 1.2;
a1 = 4;
a2 = -10;
g = 9.81;
m = 120;
L = 1.1;
R = 14;
df = 0.8;
db = 0.3;
h = 0.44;
h1 = 0.2;
b1 = 0.16;
dh = 0.33;
rd = 0.11;
rb = 0.09;
bb = 0.28;
bd = 0.0;
D = 0.028;
d = 0.6*D;
a = 0;
gamma = 0.25;
rh = 0.5*dh;



%%
% Plot for xz-plane Case 1

Nb = 908.5;
R1z = -454.25;
R2z = -454.25;

x_forces = [0, b1, L-b1, L];

F_forces = [(Nb/2), R1z, R2z, (Nb/2)];

num_forces = length(F_forces);

if length(x_forces) ~= num_forces
    error('The number of forces and positions must be the same.');
end

sum_forces = sum(F_forces);

moment_about_A = sum(F_forces .* x_forces);
R2 = moment_about_A / L;
R1 = sum_forces - R2;

x = linspace(0, L, 1000);

V = zeros(size(x));
Mz = zeros(size(x));

for i = 1:length(x)
    xi = x(i);
    V(i) = R1;
    for j = 1:num_forces
        if xi >= x_forces(j)
            V(i) = V(i) - F_forces(j);
        end
    end
    if xi >= L
        V(i) = V(i) + R2;
    end
end

% Calculate bending moment M(x)
for i = 1:length(x)
    xi = x(i);
    M_i = R1 * xi;
    for j = 1:num_forces
        if xi >= x_forces(j)
            M_i = M_i - F_forces(j) * (xi - x_forces(j));
        end
    end
    if xi >= L
        M_i = M_i + R2 * (xi - L);
    end
    Mz(i) = M_i;
end

% Plot the shear force diagram
figure;
subplot(2, 1, 1);
plot(x, V, 'b-', 'LineWidth', 2);
xlabel('Position along the beam (m)');
ylabel('Shear Force V(x) (N)');
title('Shear Force Diagram');
grid on;
xlim([0, L]);
ylim([min(V) - 10, max(V) + 10]);

% Plot the bending moment diagram
subplot(2, 1, 2);
plot(x, Mz, 'r-', 'LineWidth', 2);
xlabel('Position along the beam (m)');
ylabel('Bending Moment M(x) (Nm)');
title('Bending Moment Diagram');
grid on;
xlim([0, L]);
ylim([min(Mz) - 10, max(Mz) + 10]);

%%
% plot xy-plane Case 1

Fd = 90;
R1y = -112.5;
R2y = -112.5;
Fk = 135;

x_forces = [0, b1, L/2 ,L-b1, L];

F_forces = [(Fd/2), R1y, Fk, R2y, (Fd/2)];

num_forces = length(F_forces);

if length(x_forces) ~= num_forces
    error('The number of forces and positions must be the same.');
end

sum_forces = sum(F_forces);

moment_about_A = sum(F_forces .* x_forces);
R2 = moment_about_A / L;
R1 = sum_forces - R2;

x = linspace(0, L, 1000);

V = zeros(size(x));
My = zeros(size(x));

for i = 1:length(x)
    xi = x(i);
    V(i) = R1;
    for j = 1:num_forces
        if xi >= x_forces(j)
            V(i) = V(i) - F_forces(j);
        end
    end
    if xi >= L
        V(i) = V(i) + R2;
    end
end

% Calculate bending moment M(x)
for i = 1:length(x)
    xi = x(i);
    M_i = R1 * xi;
    for j = 1:num_forces
        if xi >= x_forces(j)
            M_i = M_i - F_forces(j) * (xi - x_forces(j));
        end
    end
    if xi >= L
        M_i = M_i + R2 * (xi - L);
    end
    My(i) = M_i;
end

% Plot the shear force diagram
figure;
subplot(2, 1, 1);
plot(x, V, 'b-', 'LineWidth', 2);
xlabel('Position along the beam (m)');
ylabel('Shear Force V(x) (N)');
title('Shear Force Diagram');
grid on;
xlim([0, L]);
ylim([min(V) - 10, max(V) + 10]);

% Plot the bending moment diagram
subplot(2, 1, 2);
plot(x, My, 'r-', 'LineWidth', 2);
xlabel('Position along the beam (m)');
ylabel('Bending Moment M(x) (Nm)');
title('Bending Moment Diagram');
grid on;
xlim([0, L]);
ylim([min(My) - 10, max(My) + 10]);

%% Nominella spänningar fall 1

x = linspace(0, L, 1000); % Creates values from 0 to L with step size 0.1

% Initialize an array for y with the same size as x
d_values = zeros(size(x));
%sigma_nom = zeros(size(x));

% Loop through each value of x and apply the conditions
for i = 1:length(x)
    if x(i) > 0 && x(i) <= b1
        d_values(i) = 0.6*D;
    elseif x(i) > b1 && x(i) <= (L-b1)
        d_values(i) = D; 
    else
        d_values(i) = 0.6*D;
    end
end

sigma_nom = zeros(size(x));

for j = 1:length(x)
    sigma_nom(j) = (32 ./ (pi * d_values(j)^3)) .* sqrt(My(j).^2 + Mz(j).^2);
end

plot(x, sigma_nom)

%%
% plot xz-plane Case 2

Nb = 1048.15;
R1z = -524;
R2z = -524;

x_forces = [0, b1, L-b1, L];

F_forces = [(Nb/2), R1z, R2z, (Nb/2)];

num_forces = length(F_forces);

if length(x_forces) ~= num_forces
    error('The number of forces and positions must be the same.');
end

sum_forces = sum(F_forces);

moment_about_A = sum(F_forces .* x_forces);
R2 = moment_about_A / L;
R1 = sum_forces - R2;

x = linspace(0, L, 1000);

V = zeros(size(x));
M = zeros(size(x));

for i = 1:length(x)
    xi = x(i);
    V(i) = R1;
    for j = 1:num_forces
        if xi >= x_forces(j)
            V(i) = V(i) - F_forces(j);
        end
    end
    if xi >= L
        V(i) = V(i) + R2;
    end
end

% Calculate bending moment M(x)
for i = 1:length(x)
    xi = x(i);
    M_i = R1 * xi;
    for j = 1:num_forces
        if xi >= x_forces(j)
            M_i = M_i - F_forces(j) * (xi - x_forces(j));
        end
    end
    if xi >= L
        M_i = M_i + R2 * (xi - L);
    end
    M(i) = M_i;
end

% Plot the shear force diagram
figure;
subplot(2, 1, 1);
plot(x, V, 'b-', 'LineWidth', 2);
xlabel('Position along the beam (m)');
ylabel('Shear Force V(x) (N)');
title('Shear Force Diagram');
grid on;
xlim([0, L]);
ylim([min(V) - 10, max(V) + 10]);

% Plot the bending moment diagram
subplot(2, 1, 2);
plot(x, M, 'r-', 'LineWidth', 2);
xlabel('Position along the beam (m)');
ylabel('Bending Moment M(x) (Nm)');
title('Bending Moment Diagram');
grid on;
xlim([0, L]);
ylim([min(M) - 10, max(M) + 10]);

%%
% plot xy-plane Case 2

Fd = 480;
R1y = -600;
R2y = -600;
Fk = 720;

x_forces = [0, b1, L/2 ,L-b1, L];

F_forces = [(Fd/2), R1y, Fk, R2y, (Fd/2)];

num_forces = length(F_forces);

if length(x_forces) ~= num_forces
    error('The number of forces and positions must be the same.');
end

sum_forces = sum(F_forces);

moment_about_A = sum(F_forces .* x_forces);
R2 = moment_about_A / L;
R1 = sum_forces - R2;

x = linspace(0, L, 1000);

V = zeros(size(x));
M = zeros(size(x));

for i = 1:length(x)
    xi = x(i);
    V(i) = R1;
    for j = 1:num_forces
        if xi >= x_forces(j)
            V(i) = V(i) - F_forces(j);
        end
    end
    if xi >= L
        V(i) = V(i) + R2;
    end
end

% Calculate bending moment M(x)
for i = 1:length(x)
    xi = x(i);
    M_i = R1 * xi;
    for j = 1:num_forces
        if xi >= x_forces(j)
            M_i = M_i - F_forces(j) * (xi - x_forces(j));
        end
    end
    if xi >= L
        M_i = M_i + R2 * (xi - L);
    end
    M(i) = M_i;
end

% Plot the shear force diagram
figure;
subplot(2, 1, 1);
plot(x, V, 'b-', 'LineWidth', 2);
xlabel('Position along the beam (m)');
ylabel('Shear Force V(x) (N)');
title('Shear Force Diagram');
grid on;
xlim([0, L]);
ylim([min(V) - 10, max(V) + 10]);

% Plot the bending moment diagram
subplot(2, 1, 2);
plot(x, M, 'r-', 'LineWidth', 2);
xlabel('Position along the beam (m)');
ylabel('Bending Moment M(x) (Nm)');
title('Bending Moment Diagram');
grid on;
xlim([0, L]);
ylim([min(M) - 10, max(M) + 10]);

%%
% plot xz-plane Case 3

Nb = 428.5;
R1z = -933.7;
R2z = -933.7;
Fb = 2035;

x_forces = [0, b1, 0.5*L-bb, 0.5*L+bb, L-b1, L];

F_forces = [(Nb/2), R1z, (1/(2*sqrt(2)))*Fb, (1/(2*sqrt(2)))*Fb, R2z, (Nb/2)];

num_forces = length(F_forces);

if length(x_forces) ~= num_forces
    error('The number of forces and positions must be the same.');
end

sum_forces = sum(F_forces);

moment_about_A = sum(F_forces .* x_forces);
R2 = moment_about_A / L;
R1 = sum_forces - R2;

x = linspace(0, L, 1000);

V = zeros(size(x));
M = zeros(size(x));

for i = 1:length(x)
    xi = x(i);
    V(i) = R1;
    for j = 1:num_forces
        if xi >= x_forces(j)
            V(i) = V(i) - F_forces(j);
        end
    end
    if xi >= L
        V(i) = V(i) + R2;
    end
end

% Calculate bending moment M(x)
for i = 1:length(x)
    xi = x(i);
    M_i = R1 * xi;
    for j = 1:num_forces
        if xi >= x_forces(j)
            M_i = M_i - F_forces(j) * (xi - x_forces(j));
        end
    end
    if xi >= L
        M_i = M_i + R2 * (xi - L);
    end
    M(i) = M_i;
end

% Plot the shear force diagram
figure;
subplot(2, 1, 1);
plot(x, V, 'b-', 'LineWidth', 2);
xlabel('Position along the beam (m)');
ylabel('Shear Force V(x) (N)');
title('Shear Force Diagram');
grid on;
xlim([0, L]);
ylim([min(V) - 10, max(V) + 10]);

% Plot the bending moment diagram
subplot(2, 1, 2);
plot(x, M, 'r-', 'LineWidth', 2);
xlabel('Position along the beam (m)');
ylabel('Bending Moment M(x) (Nm)');
title('Bending Moment Diagram');
grid on;
xlim([0, L]);
ylim([min(M) - 10, max(M) + 10]);

%%
% plot xy-plane Case 3

Fd = -1100;
R1y = 1274;
R2y = 1274;
Fb = -2035;

x_forces = [0, b1, 0.5*L-bb, 0.5*L+bb, L-b1, L];

F_forces = [(Fd/2), R1y, (1/(2*sqrt(2)))*Fb, (1/(2*sqrt(2)))*Fb , R2y, (Fd/2)];

num_forces = length(F_forces);

if length(x_forces) ~= num_forces
    error('The number of forces and positions must be the same.');
end

sum_forces = sum(F_forces);

moment_about_A = sum(F_forces .* x_forces);
R2 = moment_about_A / L;
R1 = sum_forces - R2;

x = linspace(0, L, 1000);

V = zeros(size(x));
M = zeros(size(x));

for i = 1:length(x)
    xi = x(i);
    V(i) = R1;
    for j = 1:num_forces
        if xi >= x_forces(j)
            V(i) = V(i) - F_forces(j);
        end
    end
    if xi >= L
        V(i) = V(i) + R2;
    end
end

% Calculate bending moment M(x)
for i = 1:length(x)
    xi = x(i);
    M_i = R1 * xi;
    for j = 1:num_forces
        if xi >= x_forces(j)
            M_i = M_i - F_forces(j) * (xi - x_forces(j));
        end
    end
    if xi >= L
        M_i = M_i + R2 * (xi - L);
    end
    M(i) = M_i;
end

% Plot the shear force diagram
figure;
subplot(2, 1, 1);
plot(x, V, 'b-', 'LineWidth', 2);
xlabel('Position along the beam (m)');
ylabel('Shear Force V(x) (N)');
title('Shear Force Diagram');
grid on;
xlim([0, L]);
ylim([min(V) - 10, max(V) + 10]);

% Plot the bending moment diagram
subplot(2, 1, 2);
plot(x, M, 'r-', 'LineWidth', 2);
xlabel('Position along the beam (m)');
ylabel('Bending Moment M(x) (Nm)');
title('Bending Moment Diagram');
grid on;
xlim([0, L]);
ylim([min(M) - 10, max(M) + 10]);

%%
% plot xz-plane Case 4

Vbi = 255.9;
R1z = -276.2;
R2z = -583.26;
Vby = 603.5;

x_forces = [0, b1, L-b1, L];

F_forces = [Vbi, R1z, R2z, Vby];

num_forces = length(F_forces);

if length(x_forces) ~= num_forces
    error('The number of forces and positions must be the same.');
end

sum_forces = sum(F_forces);

moment_about_A = sum(F_forces .* x_forces);
R2 = moment_about_A / L;
R1 = sum_forces - R2;

x = linspace(0, L, 1000);

V = zeros(size(x));
M = zeros(size(x));

for i = 1:length(x)
    xi = x(i);
    V(i) = R1;
    for j = 1:num_forces
        if xi >= x_forces(j)
            V(i) = V(i) - F_forces(j);
        end
    end
    if xi >= L
        V(i) = V(i) + R2;
    end
end

% Calculate bending moment M(x)
for i = 1:length(x)
    xi = x(i);
    M_i = R1 * xi;
    for j = 1:num_forces
        if xi >= x_forces(j)
            M_i = M_i - F_forces(j) * (xi - x_forces(j));
        end
    end
    if xi >= L
        M_i = M_i + R2 * (xi - L);
    end
    M(i) = M_i;
end

% Plot the shear force diagram
figure;
subplot(2, 1, 1);
plot(x, V, 'b-', 'LineWidth', 2);
xlabel('Position along the beam (m)');
ylabel('Shear Force V(x) (N)');
title('Shear Force Diagram');
grid on;
xlim([0, L]);
ylim([min(V) - 10, max(V) + 10]);

% Plot the bending moment diagram
subplot(2, 1, 2);
plot(x, M, 'r-', 'LineWidth', 2);
xlabel('Position along the beam (m)');
ylabel('Bending Moment M(x) (Nm)');
title('Bending Moment Diagram');
grid on;
xlim([0, L]);
ylim([min(M) - 10, max(M) + 10]);

%%
% plot xy-plane Case 4

Fd = 5.625;
R1y = -7.03;
Fk = 8.438;
R2y = -7.03;



x_forces = [0, b1, L/2, L-b1, L];

F_forces = [(Fd/2), R1y, Fk, R2y, (Fd/2)];

num_forces = length(F_forces);

if length(x_forces) ~= num_forces
    error('The number of forces and positions must be the same.');
end

sum_forces = sum(F_forces);

moment_about_A = sum(F_forces .* x_forces);
R2 = moment_about_A / L;
R1 = sum_forces - R2;

x = linspace(0, L, 1000);

V = zeros(size(x));
M = zeros(size(x));

for i = 1:length(x)
    xi = x(i);
    V(i) = R1;
    for j = 1:num_forces
        if xi >= x_forces(j)
            V(i) = V(i) - F_forces(j);
        end
    end
    if xi >= L
        V(i) = V(i) + R2;
    end
end

% Calculate bending moment M(x)
for i = 1:length(x)
    xi = x(i);
    M_i = R1 * xi;
    for j = 1:num_forces
        if xi >= x_forces(j)
            M_i = M_i - F_forces(j) * (xi - x_forces(j));
        end
    end
    if xi >= L
        M_i = M_i + R2 * (xi - L);
    end
    M(i) = M_i;
end

% Plot the shear force diagram
figure;
subplot(2, 1, 1);
plot(x, V, 'b-', 'LineWidth', 2);
xlabel('Position along the beam (m)');
ylabel('Shear Force V(x) (N)');
title('Shear Force Diagram');
grid on;
xlim([0, L]);
ylim([min(V) - 10, max(V) + 10]);

% Plot the bending moment diagram
subplot(2, 1, 2);
plot(x, M, 'r-', 'LineWidth', 2);
xlabel('Position along the beam (m)');
ylabel('Bending Moment M(x) (Nm)');
title('Bending Moment Diagram');
grid on;
xlim([0, L]);
ylim([min(M) - 10, max(M) + 10]);

%%
% nominell spänning case 1


I = (pi*d^4)/64;

sigma1 = 0;
sigma2 = (sqrt(75^2 + 5^2))*(b1)/I;
sigma3 = (sqrt(75^2 + 20^2))*(0.5*L)/I;
sigma4 = (sqrt(75^2 + 5^2))*(L-b1)/I;
sigma5 = 0;

x_forces = [0, b1, L/2, L-b1, L];
F_forces = [sigma1, sigma2, sigma3, sigma4, sigma5];

num_forces = length(F_forces);

if length(x_forces) ~= num_forces
    error('The number of forces and positions must be the same.');
end

sum_forces = sum(F_forces);

moment_about_A = sum(F_forces .* x_forces);
R2 = moment_about_A / L;
R1 = sum_forces - R2;

x = linspace(0, L, 1000);

V = zeros(size(x));

for i = 1:length(x)
    xi = x(i);
    V(i) = R1;
    for j = 1:num_forces
        if xi >= x_forces(j)
            V(i) = V(i) - F_forces(j);
        end
    end
    if xi >= L
        V(i) = V(i) + R2;
    end
end

% Plot the shear force diagram
figure;
plot(x, V, 'b-', 'LineWidth', 2);
xlabel('Position along the beam (m)');
ylabel('nominella spänningar (N)');
grid on;
xlim([0, L]);
ylim([min(V) - 10, max(V) + 10]);

