function fall_1(c, A, v, p, a1, a2, g, m, L, R, df, db, h, h1, b1, dh, rd, rb, bb, bd, D, d, a, gamma, rh)
    % Function to analyze shear force, bending moments, torsional moments, 
    % and nominal stresses for a beam.

    % Calculations for the xz-plane
    Nb = 908.5; % Example value
    R1z = -454.25; % Example value
    R2z = -454.25; % Example value

    x_forces_xz = [0, b1, L-b1, L];
    F_forces_xz = [(Nb/2), R1z, R2z, (Nb/2)];
    num_forces_xz = length(F_forces_xz);

    if length(x_forces_xz) ~= num_forces_xz
        error('The number of forces and positions must be the same.');
    end

    sum_forces_xz = sum(F_forces_xz);
    moment_about_A_xz = sum(F_forces_xz .* x_forces_xz);
    R2_xz = moment_about_A_xz / L;
    R1_xz = sum_forces_xz - R2_xz;

    x = linspace(0, L, 1000);
    V_xz = zeros(size(x));
    Mz = zeros(size(x));

    for i = 1:length(x)
        xi = x(i);
        V_xz(i) = R1_xz;
        for j = 1:num_forces_xz
            if xi >= x_forces_xz(j)
                V_xz(i) = V_xz(i) - F_forces_xz(j);
            end
        end
        if xi >= L
            V_xz(i) = V_xz(i) + R2_xz;
        end
    end

    for i = 1:length(x)
        xi = x(i);
        M_i = R1_xz * xi;
        for j = 1:num_forces_xz
            if xi >= x_forces_xz(j)
                M_i = M_i - F_forces_xz(j) * (xi - x_forces_xz(j));
            end
        end
        if xi >= L
            M_i = M_i + R2_xz * (xi - L);
        end
        Mz(i) = M_i;
    end

    % Plot shear force and bending moment for xz-plane
    figure;
    subplot(2, 1, 1);
    plot(x, V_xz, 'b-', 'LineWidth', 2);
    xlabel('Position along the beam (m)');
    ylabel('Shear Force V(x) (N)');
    title('Shear Force Diagram (xz-plane)');
    grid on;
    xlim([0, L]);
    ylim([min(V_xz) - 10, max(V_xz) + 10]);

    subplot(2, 1, 2);
    plot(x, Mz, 'r-', 'LineWidth', 2);
    xlabel('Position along the beam (m)');
    ylabel('Bending Moment M(x) (Nm)');
    title('Bending Moment Diagram (xz-plane)');
    grid on;
    xlim([0, L]);
    ylim([min(Mz) - 10, max(Mz) + 10]);

    % Calculations for the xy-plane
    Fd = 90;
    R1y = -112.5;
    R2y = -112.5;
    Fk = 135;

    x_forces_xy = [0, b1, L/2, L-b1, L];
    F_forces_xy = [(Fd/2), R1y, Fk, R2y, (Fd/2)];
    num_forces_xy = length(F_forces_xy);

    if length(x_forces_xy) ~= num_forces_xy
        error('The number of forces and positions must be the same.');
    end

    sum_forces_xy = sum(F_forces_xy);
    moment_about_A_xy = sum(F_forces_xy .* x_forces_xy);
    R2_xy = moment_about_A_xy / L;
    R1_xy = sum_forces_xy - R2_xy;

    V_xy = zeros(size(x));
    My = zeros(size(x));

    for i = 1:length(x)
        xi = x(i);
        V_xy(i) = R1_xy;
        for j = 1:num_forces_xy
            if xi >= x_forces_xy(j)
                V_xy(i) = V_xy(i) - F_forces_xy(j);
            end
        end
        if xi >= L
            V_xy(i) = V_xy(i) + R2_xy;
        end
    end

    for i = 1:length(x)
        xi = x(i);
        M_i = R1_xy * xi;
        for j = 1:num_forces_xy
            if xi >= x_forces_xy(j)
                M_i = M_i - F_forces_xy(j) * (xi - x_forces_xy(j));
            end
        end
        if xi >= L
            M_i = M_i + R2_xy * (xi - L);
        end
        My(i) = M_i;
    end

    % Plot shear force and bending moment for xy-plane
    figure;
    subplot(2, 1, 1);
    plot(x, V_xy, 'b-', 'LineWidth', 2);
    xlabel('Position along the beam (m)');
    ylabel('Shear Force V(x) (N)');
    title('Shear Force Diagram (xy-plane)');
    grid on;
    xlim([0, L]);
    ylim([min(V_xy) - 10, max(V_xy) + 10]);

    subplot(2, 1, 2);
    plot(x, My, 'r-', 'LineWidth', 2);
    xlabel('Position along the beam (m)');
    ylabel('Bending Moment M(x) (Nm)');
    title('Bending Moment Diagram (xy-plane)');
    grid on;
    xlim([0, L]);
    ylim([min(My) - 10, max(My) + 10]);

    % Torsional moment calculations
    Moment = [0.5*Fd*rh, -Fk*rd+0.5*Fd*rh];
    moment_verkan = [0, L/2];

    num_moment = length(Moment);
    sum_moment = sum(Moment);

    moment_about_A = sum(Moment .* moment_verkan);
    R2 = moment_about_A / L;
    R1 = sum_moment - R2;

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

    % Plot torsional moment
    figure;
    plot(x, V, 'b-', 'LineWidth', 2);
    xlabel('Position along the beam (m)');
    ylabel('Torsional Moment (N)');
    title("Torsional Moment");
    grid on;
    xlim([0, L]);
    ylim([min(V) - 10, max(V) + 10]);

    % Nominal stresses
    d_values = zeros(size(x));
    sigma_nom = zeros(size(x));

    for i = 1:length(x)
        if x(i) > 0 && x(i) <= b1
            d_values(i) = 0.6 * D;
        elseif x(i) > b1 && x(i) <= (L - b1)
            d_values(i) = D;
        else
            d_values(i) = 0.6 * D;
        end
    end

    for j = 1:length(x)
        sigma_nom(j) = (32 ./ (pi * d_values(j).^3)) .* sqrt(My(j).^2 + Mz(j).^2);
    end

    % Plot nominal stresses
    figure;
    plot(x, sigma_nom, 'k-', 'LineWidth', 2);
    xlabel('Position along the beam (m)');
    ylabel('Nominal Stress \sigma_{nom} (Pa)');
    title('Nominal Stress Distribution');
    grid on;
    xlim([0, L]);

    fall_1(0.3, 0.45, 120/3.6, 1.2, 4, -10, 9.81, 120, 1.1, 14, 0.8, 0.3, 0.44, 0.2, 0.16, 0.33, 0.11, 0.09, 0.28, 0, 0.028, 0.6*0.028, 0, 0.25, 0.5*0.33)
end
