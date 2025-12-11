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
% Vridmoment Case 1

Fd = 90;
Fk = 135;
rh = 0.5*0.33;
rd = 0.11;
L = 1.1;

Moment = [0.5*Fd*rh, -Fk*rd+0.5*Fd*rh];
moment_verkan = [0, L/2];

num_moment = length(Moment);
sum_moment = sum(Moment);

moment_about_A = sum(Moment .* moment_verkan);
R2 = moment_about_A / L;
R1 = sum_moment - R2;

x = linspace(0, L, 1000);

V = zeros(size(x));

for i = 1:length(x)
    xi = x(i);
    V(i) = R1;
    for j = 1:num_moment
        if xi >= moment_verkan(j)
            V(i) = V(i) - Moment(j);
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
ylabel('Vridande moment (N)');
title("Fall 1: Vridande Moment")
grid on;
xlim([0, L]);
ylim([min(V) - 10, max(V) + 10]);

%%
% Vridmoment Case 2 

Fd = 480;
Fk = 720;
rh = 0.5*0.33;
rd = 0.11;
L = 1.1;

Moment = [Fd*rh/2, -Fk*rd+0.5*Fd*rh];
moment_verkan = [0, L/2];

num_moment = length(Moment);
sum_moment = sum(Moment);

moment_about_A = sum(Moment .* moment_verkan);
R2 = moment_about_A / L;
R1 = sum_moment - R2;

x = linspace(0, L, 1000);

V = zeros(size(x));

for i = 1:length(x)
    xi = x(i);
    V(i) = R1;
    for j = 1:num_moment
        if xi >= moment_verkan(j)
            V(i) = V(i) - Moment(j);
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
ylabel('Vridande moment (N)');
title("Fall 2: Vridande Moment")
grid on;
xlim([0, L]);
ylim([min(V) - 10, max(V) + 10]);

%%
% Vridmoment fall 3

Fd = -1100;
Fb = -2035;
rh = 0.5*0.33;
rb = 0.09;
L = 1.1;

Moment = [-0.5*Fd*rh, (-0.5*Fd*rh) - (Fb*rb), -0.5*Fd*rh - 2*Fb*rb];
moment_verkan = [0, 0.5*L-bb, 0.5*L+bb];

num_moment = length(Moment);
sum_moment = sum(Moment);

moment_about_A = sum(Moment .* moment_verkan);
R2 = moment_about_A / L;
R1 = sum_moment - R2;

x = linspace(0, L, 1000);

V = zeros(size(x));

for i = 1:length(x)
    xi = x(i);
    V(i) = R1;
    for j = 1:num_moment
        if xi >= moment_verkan(j)
            V(i) = V(i) - Moment(j);
        end
    end
    if xi >= L
        V(i) = V(i) + R2;
    end
end

% Plot the shear force diagram
figure;
plot(x, V, 'b-', 'LineWidth', 2);
xlabel('Position (m)');
ylabel('Vridande moment (N)');
title("Fall 3 Vridande moment")
grid on;
xlim([0, L]);
ylim([min(V) - 20, max(V) + 20]);

%%
% Normalspänning fall 4

Hbi = 128.89;
Rx = 432.9;
Hby = 304;

Moment = [-Rx-Hby ,Hby];
moment_verkan = [0, b1];

num_moment = length(Moment);
sum_moment = sum(Moment);

moment_about_A = sum(Moment .* moment_verkan);
R2 = moment_about_A / L;
R1 = sum_moment - R2;

x = linspace(0, L, 1000);

V = zeros(size(x));

for i = 1:length(x)
    xi = x(i);
    V(i) = R1;
    for j = 1:num_moment
        if xi >= moment_verkan(j)
            V(i) = V(i) - Moment(j);
        end
    end
    if xi >= L
        V(i) = V(i) + R2;
    end
end

% Plot the shear force diagram
figure;
plot(x, V, 'b-', 'LineWidth', 2);
xlabel('Position (m)');
ylabel('Normalkraft (N)');
title("Fall 4: Normalkraft")
grid on;
xlim([0, L]);
ylim([min(V) - 20, max(V) + 20]);

%% Vridande moment 4

Fd = 5.625;
Fk = 8.4375;
rh = 0.5*0.33;
rd = 0.11;
L = 1.1;

Moment = [-Fd*rh/2, Fk*rd-0.5*Fd*rh];
moment_verkan = [0, L/2];

num_moment = length(Moment);
sum_moment = sum(Moment);

moment_about_A = sum(Moment .* moment_verkan);
R2 = moment_about_A / L;
R1 = sum_moment - R2;

x = linspace(0, L, 1000);

V = zeros(size(x));

for i = 1:length(x)
    xi = x(i);
    V(i) = R1;
    for j = 1:num_moment
        if xi >= moment_verkan(j)
            V(i) = V(i) - Moment(j);
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
ylabel('Vridande moment (N)');
grid on;
xlim([0, L]);
ylim([min(V), max(V)]);