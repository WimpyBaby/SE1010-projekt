clear; clc; close all;

% Relevant Parameters
L = 1.1;        % Beam length
b1 = 0.16;      % Position parameter
D = 0.055;      % Diameter used for stress calculations
dh = 0.33;      % Used to compute rh
rh = 0.5 * dh;  % Half of dh
rd = 0.11;      % Used in Mv calculations
bb = 0.28;      % Used for certain force positions
K_tb = 1.7;
K_tv = 1.5;
sigma_lim = 110*10^6;

%--------------------------------------------------------------------------
% Fall 4: XZ-Plane
%--------------------------------------------------------------------------

% Given forces
Vbi = 255.9;
R1z = -276.2;
R2z = -583.26;
Vby = 603.5;

x_forces_xz = [0, b1, L - b1, L];
F_forces_xz = [Vbi, R1z, R2z, Vby];

if length(x_forces_xz) ~= length(F_forces_xz)
    error('Number of XZ forces must match number of positions.');
end

sum_forces_xz = sum(F_forces_xz);
moment_about_A_xz = sum(F_forces_xz .* x_forces_xz);
R2_xz = moment_about_A_xz / L;
R1_xz = sum_forces_xz - R2_xz;

x = linspace(0, L, 1000);
V_xz = zeros(size(x));  % Shear force (XZ)
Mz = zeros(size(x));    % Bending moment (XZ)

% Calculate shear force V(x) in XZ-plane
for i = 1:length(x)
    xi = x(i);
    V_xz(i) = R1_xz;
    for j = 1:length(F_forces_xz)
        if xi >= x_forces_xz(j)
            V_xz(i) = V_xz(i) - F_forces_xz(j);
        end
    end
    if xi >= L
        V_xz(i) = V_xz(i) + R2_xz;
    end
end

% Calculate bending moment Mz(x) in XZ-plane
for i = 1:length(x)
    xi = x(i);
    Mz_i = R1_xz * xi;
    for j = 1:length(F_forces_xz)
        if xi >= x_forces_xz(j)
            Mz_i = Mz_i - F_forces_xz(j) * (xi - x_forces_xz(j));
        end
    end
    if xi >= L
        Mz_i = Mz_i + R2_xz * (xi - L);
    end
    Mz(i) = Mz_i;
end

%--------------------------------------------------------------------------
% Fall 4: XY-Plane
%--------------------------------------------------------------------------

Fd = 5.625;
R1y = -7.03;
Fk = 8.438;
R2y = -7.03;

x_forces_xy = [0, b1, L/2, L - b1, L];
F_forces_xy = [(Fd/2), R1y, Fk, R2y, (Fd/2)];

if length(x_forces_xy) ~= length(F_forces_xy)
    error('Number of XY forces must match number of positions.');
end

sum_forces_xy = sum(F_forces_xy);
moment_about_A_xy = sum(F_forces_xy .* x_forces_xy);
R2_xy = moment_about_A_xy / L;
R1_xy = sum_forces_xy - R2_xy;

V_xy = zeros(size(x)); % Shear force (XY)
My = zeros(size(x));   % Bending moment (XY)

% Calculate shear force V(x) in XY-plane
for i = 1:length(x)
    xi = x(i);
    V_xy(i) = R1_xy;
    for j = 1:length(F_forces_xy)
        if xi >= x_forces_xy(j)
            V_xy(i) = V_xy(i) - F_forces_xy(j);
        end
    end
    if xi >= L
        V_xy(i) = V_xy(i) + R2_xy;
    end
end

% Calculate bending moment My(x) in XY-plane
for i = 1:length(x)
    xi = x(i);
    My_i = R1_xy * xi;
    for j = 1:length(F_forces_xy)
        if xi >= x_forces_xy(j)
            My_i = My_i - F_forces_xy(j) * (xi - x_forces_xy(j));
        end
    end
    if xi >= L
        My_i = My_i + R2_xy * (xi - L);
    end
    My(i) = My_i;
end

%--------------------------------------------------------------------------
% Fall 4: Normalkraft
%--------------------------------------------------------------------------

Hbi = 128.89;
Rx = 432.9;
Hby = 304;

Moment_norm = [-Rx - Hby, Hby];
moment_verkan_norm = [0, b1];

if length(Moment_norm) ~= length(moment_verkan_norm)
    error('Number of normal moments must match positions.');
end

sum_moment_norm = sum(Moment_norm);
moment_about_A_norm = sum(Moment_norm .* moment_verkan_norm);
R2_norm = moment_about_A_norm / L;
R1_norm = sum_moment_norm - R2_norm;

Norm = zeros(size(x));  % Normal force distribution

for i = 1:length(x)
    xi = x(i);
    Norm(i) = R1_norm;
    for j = 1:length(Moment_norm)
        if xi >= moment_verkan_norm(j)
            Norm(i) = Norm(i) - Moment_norm(j);
        end
    end
    if xi >= L
        Norm(i) = Norm(i) + R2_norm;
    end
end

%--------------------------------------------------------------------------
% Fall 4: Nominella Spänningar
%--------------------------------------------------------------------------

d_values = zeros(size(x));
for i = 1:length(x)
    if x(i) > 0 && x(i) <= b1
        d_values(i) = 0.6 * D;
    elseif x(i) > b1 && x(i) <= (L - b1)
        d_values(i) = D; 
    else
        d_values(i) = 0.6 * D;
    end
end

sigma_nom = zeros(size(x));
for j = 1:length(x)
    sigma_nom(j) = (4 * Norm / (pi * d_values.^2)) + (32 / (pi * d_values(j)^3)) * sqrt(My(j)^2 + Mz(j)^2);
end

%--------------------------------------------------------------------------
% Fall 4: Vridande Moment (Mv)
%--------------------------------------------------------------------------

Fd = 5.625;   % Reconfirming Fd for this section
Fk = 8.4375;  % Slightly different precision from above
rh = 0.5 * 0.33; % Redefining rh for clarity, same as before
rd = 0.11;    % Already defined above, kept for clarity
L = 1.1;      % Reconfirm L

Moment_mv = [-Fd*rh/2, Fk*rd - 0.5*Fd*rh];
moment_verkan_mv = [0, L/2];

if length(Moment_mv) ~= length(moment_verkan_mv)
    error('Number of torsional moments must match positions.');
end

sum_moment_mv = sum(Moment_mv);
moment_about_A_mv = sum(Moment_mv .* moment_verkan_mv);
R2_mv = moment_about_A_mv / L;
R1_mv = sum_moment_mv - R2_mv;

Mv = zeros(size(x));
for i = 1:length(x)
    xi = x(i);
    Mv(i) = R1_mv;
    for j = 1:length(Moment_mv)
        if xi >= moment_verkan_mv(j)
            Mv(i) = Mv(i) - Moment_mv(j);
        end
    end
    if xi >= L
        Mv(i) = Mv(i) + R2_mv;
    end
end

%--------------------------------------------------------------------------
% Fall 4: Skjuvspänning (tau_v) due to Mv
%--------------------------------------------------------------------------

tau_v = zeros(size(x));
for i = 1:length(x)
    tau_v(i) = (32 / (pi * d_values(i)^3)) * Mv(i);
end

tau_y = zeros(size(x));
tau_z = zeros(size(x));
for i = 1:length(x)
    tau_y(i) = 16 .* My(i) ./ (pi .* d_values(i).^3);
    tau_z(i) = 16 .* Mz(i) ./ (pi .* d_values(i).^3);
end

% -------------------------------------------------------------------------
% effektivspänning (sigma_e)
% -------------------------------------------------------------------------
sigma_e = zeros(size(x));
for i = 1:length(x)
    sigma_e(i) = sqrt(sigma_nom(i)^2 + 3 * tau_y(i)^2 + 3 * tau_z(i)^2);
end

%--------------------------------------------------------------------------
% sigma max (sigma_max)
%--------------------------------------------------------------------------

sigma_max = zeros(size(x));
for i = 1:length(x)
    sigma_max(i) = sqrt(K_tb * sigma_nom(i)^2 + 3 * (K_tv * tau_y(i)^2 + K_tv * tau_z(i)^2));
end

%--------------------------------------------------------------------------
% Combine All Plots into One Figure
%
% Using a 4x2 grid for the 8 plots:
% 1: XZ Shear Force (V_xz)
% 2: XZ Bending Moment (Mz)
% 3: XY Shear Force (V_xy)
% 4: XY Bending Moment (My)
% 5: Normalkraft (Norm)
% 6: Nominella Spänningar (sigma_nom)
% 7: Vridande Moment (Mv)
% 8: Skjuvspänning (tau_v)
%--------------------------------------------------------------------------

figure;

sgtitle('Fall 4');

subplot(5,2,1);
plot(x, V_xz, 'b-', 'LineWidth', 2);
xlabel('Position (m)');
ylabel('N');
xlim([0, L]);
ylim([min(V_xz), max(V_xz)]);
title('Fall 4: XZ Shear Force');
grid on;

subplot(5,2,2);
plot(x, Mz, 'r-', 'LineWidth', 2);
xlabel('Position (m)');
ylabel('Nm');
xlim([0, L]);
ylim([min(Mz), max(Mz)]);
title('Fall 4: XZ Bending Moment');
grid on;

subplot(5,2,3);
plot(x, V_xy, 'b-', 'LineWidth', 2);
xlabel('Position (m)');
ylabel('N');
xlim([0, L]);
ylim([min(V_xy), max(V_xy)]);
title('Fall 4: XY Shear Force');
grid on;

subplot(5,2,4);
plot(x, My, 'r-', 'LineWidth', 2);
xlabel('Position (m)');
ylabel('Nm');
xlim([0, L]);
ylim([min(My), max(My)]);
title('Fall 4: XY Bending Moment');
grid on;

subplot(5,2,5);
plot(x, Norm, 'k-', 'LineWidth', 2);
xlabel('Position (m)');
ylabel('N');
xlim([0, L]);
ylim([min(Norm), max(Norm)]);
title('Normalkraft');
grid on;

subplot(5,2,6);
plot(x, sigma_nom, 'b-', 'LineWidth', 2);
xlabel('Position (m)');
ylabel('Pa');
xlim([0, L]);
ylim([min(sigma_nom), max(sigma_nom)]);
title('\sigma_{nom}');
grid on;

subplot(5,2,7);
plot(x, Mv, 'r-', 'LineWidth', 2);
xlabel('Position (m)');
ylabel('Nm');
xlim([0, L]);
ylim([min(Mv), max(Mv)]);
title('M_v');
grid on;

subplot(5,2,8);
plot(x, tau_v, 'k-', 'LineWidth', 2);
ylabel('Pa');
xlim([0, L]);
ylim([min(tau_v), max(tau_v)]);
title('\tau_v');
grid on;

subplot(5,2,9);
plot(x, sigma_e, 'k-', 'LineWidth', 2);
ylabel('Pa');
xlim([0, L]);
ylim([min(sigma_e), max(sigma_e)]);
title('\sigma_e');
grid on;

subplot(5,2,10);
plot(x, sigma_max, 'k-', 'LineWidth', 2);
hold on;
yline(sigma_lim, 'g', 'LineWidth', 2);
xlabel('Position (m)');
ylabel('Pa');
xlim([0, L]);

allSigmaValues = [sigma_max(:); sigma_lim];
ylim([min(allSigmaValues), max(allSigmaValues) + 1e7]);
title('\sigma_{max}');
grid on;

legend('\sigma_{max}', '\sigma_{lim}');

