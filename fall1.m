clear; clc; close all;

% Common Parameters
L = 1.1;       % Beam length
b1 = 0.16;     % Position parameter along the beam
D = 0.028;     % Diameter for nominal stress calculations
K_tb = 1.7;
K_tv = 1.5;


%--------------------------------------------------------------------------
% Fall 1: XZ-Plane
%--------------------------------------------------------------------------

Nb = 908.5;
R1z = -454.25;
R2z = -454.25;

x_forces_xz = [0, b1, L - b1, L];
F_forces_xz = [(Nb/2), R1z, R2z, (Nb/2)];

if length(x_forces_xz) ~= length(F_forces_xz)
    error('Number of forces must equal number of positions in XZ scenario.');
end

sum_forces_xz = sum(F_forces_xz);
moment_about_A_xz = sum(F_forces_xz .* x_forces_xz);
R2_xz = moment_about_A_xz / L;
R1_xz = sum_forces_xz - R2_xz;

x = linspace(0, L, 1000);
V_xz = zeros(size(x));  % Shear force distribution in XZ
Mz = zeros(size(x));    % Bending moment distribution in XZ

% Calculate V(x) for XZ
for i = 1:length(x)
    xi = x(i);
    V_xz(i) = R1_xz;
    for j = 1:length(F_forces_xz)
        if xi >= x_forces_xz(j)
            V_xz(i) = V_xz(i) - F_forces_xz(j);
        end
    end
    if xi >= L
        V_xz(i) = V_xz(i) + R2_xz;sigma_lim = 110*10^6;
    end
end

% Calculate Mz(x) for XZ
for i = 1:length(x)
    xi = x(i);
    M_i = R1_xz * xi;
    for j = 1:length(F_forces_xz)
        if xi >= x_forces_xz(j)
            M_i = M_i - F_forces_xz(j) * (xi - x_forces_xz(j));
        end
    end
    if xi >= L
        M_i = M_i + R2_xz * (xi - L);
    end
    Mz(i) = M_i;
end

%--------------------------------------------------------------------------
% Fall 1: XY-Plane
%--------------------------------------------------------------------------

Fd = 90;
R1y = -112.5;
R2y = -112.5;
Fk = 135;

x_forces_xy = [0, b1, L/2, L - b1, L];
F_forces_xy = [(Fd/2), R1y, Fk, R2y, (Fd/2)];

if length(x_forces_xy) ~= length(F_forces_xy)
    error('Number of forces must equal number of positions in XY scenario.');
end

sum_forces_xy = sum(F_forces_xy);
moment_about_A_xy = sum(F_forces_xy .* x_forces_xy);
R2_xy = moment_about_A_xy / L;
R1_xy = sum_forces_xy - R2_xy;

V_xy = zeros(size(x));  % Shear force distribution in XY
My = zeros(size(x));    % Bending moment distribution in XY

% Calculate V(x) for XY
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

% Calculate My(x) for XY
for i = 1:length(x)
    xi = x(i);
    M_i = R1_xy * xi;
    for j = 1:length(F_forces_xy)
        if xi >= x_forces_xy(j)
            M_i = M_i - F_forces_xy(j) * (xi - x_forces_xy(j));
        end
    end
    if xi >= L
        M_i = M_i + R2_xy * (xi - L);
    end
    My(i) = M_i;
end

%--------------------------------------------------------------------------
% Nominella Spänningar
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
    sigma_nom(j) = (32 / (pi * d_values(j)^3)) * sqrt(My(j)^2 + Mz(j)^2);
end

%--------------------------------------------------------------------------
% Vridande Moment (Mv)
%--------------------------------------------------------------------------

rh = 0.5 * 0.33;
rd = 0.11;

Moment = [0.5*Fd*rh, -Fk*rd + 0.5*Fd*rh];
moment_verkan = [0, L/2];

if length(Moment) ~= length(moment_verkan)
    error('Number of moments must match number of moment positions.');
end

sum_moment = sum(Moment);
moment_about_A = sum(Moment .* moment_verkan);
R2_v = moment_about_A / L;
R1_v = sum_moment - R2_v;

Mv = zeros(size(x));
for i = 1:length(x)
    xi = x(i);
    Mv(i) = R1_v;
    for j = 1:length(Moment)
        if xi >= moment_verkan(j)
            Mv(i) = Mv(i) - Moment(j);
        end
    end
    if xi >= L
        Mv(i) = Mv(i) + R2_v;
    end
end

%--------------------------------------------------------------------------
% Skjuvspänning (tau_v)
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

%--------------------------------------------------------------------------
% effektivspänning (sigma_e)
%--------------------------------------------------------------------------

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

% -------------------------------------------------------------------------
% Combine all plots into one figure
% We will use a 4x2 grid of subplots:
%
% 1: XZ Shear Force (V_xz)
% 2: XZ Bending Moment (Mz)
% 3: XY Shear Force (V_xy)
% 4: XY Bending Moment (My)
% 5: Nominal Stress (sigma_nom)
% 6: Torsional Moment (Mv)
% 7: Shear Stress (tau_v)
% 8: Unused
% -------------------------------------------------------------------------
sgtitle('Fall 1')

subplot(5,2,1);
plot(x, V_xz, 'b-', 'LineWidth', 2);
ylabel('N');
xlim([0, L]);
ylim([min(V_xz), max(V_xz)]);
title('Fall 2: XZ Shear Force');
grid on;

subplot(5,2,2);
plot(x, Mz, 'r-', 'LineWidth', 2);
ylabel('Nm');
xlim([0, L]);
ylim([min(Mz), max(Mz)]);
title('Fall 2: XZ Bending Moment');
grid on;

subplot(5,2,3);
plot(x, V_xy, 'b-', 'LineWidth', 2);
ylabel('N');
xlim([0, L]);
ylim([min(V_xy), max(V_xy)]);
title('Fall 2: XY Shear Force');
grid on;

subplot(5,2,4);
plot(x, My, 'r-', 'LineWidth', 2);
ylabel('Nm');
xlim([0, L]);
ylim([min(My), max(My)]);
title('Fall 2: XY Bending Moment');
grid on;

subplot(5,2,5);
plot(x, sigma_nom, 'b-', 'LineWidth', 2);
ylabel('Pa');
xlim([0, L]);
ylim([min(sigma_nom), max(sigma_nom)]);
title('Fall 2: Nominella Spänningar');
grid on;

subplot(5,2,6);
plot(x, Mv, 'r-', 'LineWidth', 2);
ylabel('Nm');
xlim([0, L]);
ylim([min(Mv), max(Mv)]);
title('M_v');
grid on;

subplot(5,2,7); 
plot(x, tau_v, 'k-', 'LineWidth', 2);
ylabel('Pa');
xlim([0, L]);
ylim([min(tau_v), max(tau_v)]);
title('\tau_v');
grid on;

subplot(5,2,8)
plot(x, sigma_e, 'k-', 'LineWidth', 2);
ylabel('Pa');
xlim([0, L]);
ylim([min(sigma_e), max(sigma_e)]);
title('\sigma_e');
grid on;

subplot(5,2,[9,10]);
plot(x, sigma_max, 'k-', 'LineWidth', 2);
hold on
yline(sigma_lim, 'g', 'LineWidth', 2);
xlabel('Position (m)');
ylabel('Pa');
title('\sigma_{max}');
xlim([0, L]);

allSigmaValues = [sigma_max(:); sigma_lim];
ylim([min(allSigmaValues), max(allSigmaValues) + 1e7]);

title('\sigma_{max}');
grid on;
legend('\sigma_{max}', '\sigma_{lim}');
