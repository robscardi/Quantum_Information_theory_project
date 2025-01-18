clear variables
close all
g = zeros(9, 30);
q = zeros(9, 30);

x = 2:10;

for i=2:10
a = load("Figures\"+ i +"_bit_0_disp_9_km.mat");
g(i,:) = a(1).MI_vector;
a = load("Figures\QBER_"+ i +"_bit_0_disp_9_km.mat");
q(i,:) = a(1).QBER;
end

for j = 1:length(g(:,1))
    y(j) = mean(g(j, 30-4:end));
end
y = y(2:end);
% Plot the data
figure;
plot(x, y, 'o-', 'LineWidth', 2);
hold on;

% Loop through consecutive points to highlight y-axis differences
for i = 1:length(x)-1
    x1 = x(i);
    y1 = y(i);
    x2 = x(i+1);
    y2 = y(i+1);
    
    % Plot a vertical dashed line between each pair of consecutive points
    plot([x1, x1], [y1, y2], 'r--', 'LineWidth', 2);
    
    % Annotate the difference in y-values
    dy = abs(y2 - y1);
    text(x1 + 0.06, (y1 + y2) / 2 +0.1, ['\Delta y = ' num2str(dy)], 'Color', 'r', 'FontSize', 12);
end

% Optionally, add labels
xlabel('X-axis');
ylabel('Y-axis');
title('Highlight Y-Axis Differences Between Consecutive Points');
grid on;

hold off;

