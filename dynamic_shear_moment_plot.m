function dynamic_shear_moment_plot_5_forces()

    % Create a UI figure for the slider and the plot
    fig = uifigure('Name', 'Dynamic Shear and Moment Plot (5 Forces)');
    
    % Add a slider to control the position of the fourth force (a)
    sld = uislider(fig, 'Value', 2, 'Limits', [0, 10], 'Position', [100, 50, 400, 3], ...
                   'ValueChanging', @(sld,event) sliderMoving(event, fig));
    
    % Initial plot setup
    ax = uiaxes(fig, 'Position', [100, 150, 400, 300]);
    ax.XLabel.String = 'Position along the axle (m)';
    ax.YLabel.String = 'Force / Moment';
    title(ax, 'Shear Force and Bending Moment Diagrams');
    
    % Call the function once at initialization with the initial slider value
    sliderMoving(struct('Value', sld.Value), fig);
    
    % Callback function to update the plot based on slider position
    function sliderMoving(event, fig)
        % Axle parameters
        L = 10;                % Length of the axle (meters)
        
        % Forces on the axle
        F1 = 500;              % Force at the far left (N)
        F2 = 700;              % Force at the far right (N)
        F3 = 900;              % Force at the fixed position b1 (N)
        F4 = 1000;             % Movable force (N) controlled by slider
        F5 = 600;              % Another force at a fixed position x5
        
        % Positions of the forces
        b1 = 6;                % Fixed position of F3 (meters)
        a = event.Value;       % Movable force's position (controlled by slider)
        x5 = 8;                % Position of F5 (meters)
        
        % Reactions at supports (assume axle is simply supported)
        % Static equilibrium: sum of forces = 0, sum of moments = 0
        
        % Calculate reactions at supports (left and right ends)
        % Use moments about one of the supports and force equilibrium to get R1 and R2
        % Summing moments about the left support:
        R1 = (F2*L + F3*(b1) + F4*a + F5*x5) / L;  % Reaction at left support (x = 0)
        R2 = F1 + F2 + F3 + F4 + F5 - R1;          % Reaction at right support (x = L)
        
        % Discretize the axle
        x = linspace(0, L, 1000);
        
        % Shear Force Calculation
        V = zeros(size(x));
        V(x <= 0) = R1;                 % Left support to far left
        V(x > 0 & x <= a) = R1 - F1;    % Shear from far left to movable force at 'a'
        V(x > a & x <= b1) = R1 - F1 - F4;  % Shear from movable force 'a' to fixed position b1
        V(x > b1 & x <= x5) = R1 - F1 - F4 - F3; % Shear from fixed point to F5
        V(x > x5) = R1 - F1 - F4 - F3 - F5;   % Shear after F5
        
        % Bending Moment Calculation
        M = zeros(size(x));
        M(x <= 0) = 0;                               % Left support
        M(x > 0 & x <= a) = R1 .* x(x > 0 & x <= a); % Bending moment up to a
        M(x > a & x <= b1) = R1 .* x(x > a & x <= b1) - F1 .* (x(x > a & x <= b1));  % Moment after movable force
        M(x > b1 & x <= x5) = R1 .* x(x > b1 & x <= x5) - F1 .* (x(x > b1 & x <= x5)) - F4 .* (x(x > b1 & x <= x5) - a);
        M(x > x5) = M(x > b1 & x <= x5) - F5 .* (x(x > x5) - x5); % Moment after F5
        
        % Plot the Shear Force and Moment on the same axes
        ax = findobj(fig, 'Type', 'axes');
        cla(ax);  % Clear previous plot
        
        yyaxis(ax, 'left');
        plot(ax, x, V, 'b', 'LineWidth', 2);
        ylabel(ax, 'Shear Force (N)');
        
        yyaxis(ax, 'right');
        plot(ax, x, M, 'r', 'LineWidth', 2);
        ylabel(ax, 'Bending Moment (Nm)');
        
        % Update plot titles and labels
        ax.XLabel.String = 'Position along the axle (m)';
        title(ax, ['Shear and Moment Diagrams (Force at Position: ', num2str(a), ' m)']);
        
    end

end
