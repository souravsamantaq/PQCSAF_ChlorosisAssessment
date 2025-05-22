function display_score(user_value)
    % Validate user input
    if user_value < 1 || user_value > 4
        error('Value must be between 1 and 4');
    end

    % Parameters
    min_val = 1;
    max_val = 4;
    n_ticks = max_val - min_val + 1;

    % Gauge configuration
    theta = linspace(pi, 0, 100); % half-circle
    r = 1;
    x = r * cos(theta);
    y = r * sin(theta);

    % Plot gauge arc
    figure;
    plot(x, y, 'k', 'LineWidth', 4);
    hold on;
    axis equal;
    axis off;

    % Tick marks and labels
    for i = 0:n_ticks-1
        val = min_val + i;
        angle = pi - (i / (n_ticks - 1)) * pi;
        x_tick = 1.1 * cos(angle);
        y_tick = 1.1 * sin(angle);
        text(x_tick, y_tick, num2str(val), ...
            'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'middle', ...
            'FontSize', 12, 'FontWeight', 'bold');
    end

    % Draw needle
    % Map user_value to angle
    angle = pi - ((user_value - min_val) / (max_val - min_val)) * pi;
    needle_length = 0.8;
    x_needle = [0, needle_length * cos(angle)];
    y_needle = [0, needle_length * sin(angle)];
    %figure(2)
    plot(x_needle, y_needle, 'r', 'LineWidth', 3);
    plot(0, 0, 'ko', 'MarkerFaceColor', 'k'); % needle pivot

    % Title
    title(sprintf('Severity score: %.2f', user_value), 'FontSize', 14);
end
