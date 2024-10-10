g = zeros(6, length(x));
q = zeros(6, length(x));

for i=2:6
a = load("Figures\"+ i +"_bit_0_disp_9_km.mat");
g(i,:) = a(1).MI_vector;
a = load("Figures\QBER_"+ i +"_bit_0_disp_9_km.mat");
q(i,:) = a(1).QBER;
end

C_h = ((max_phot+1).*log2(max_phot+1)- (max_phot).*log2(max_phot))*0.5;
C_s2 = log2(1+2*max_phot)*0.5;

for j = 1:length(g(:,1))
    x(j) = mean(g(j, 30-4:end));
end
y = 2:10;
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
    text(x1 + 0.2, (y1 + y2) / 2, ['\Delta y = ' num2str(dy)], 'Color', 'r', 'FontSize', 12);
end

% Optionally, add labels
xlabel('X-axis');
ylabel('Y-axis');
title('Highlight Y-Axis Differences Between Consecutive Points');
grid on;

hold off;

