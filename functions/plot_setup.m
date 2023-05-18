function plot_setup(scattererPoints, txPoints, rxPoints)
%PLOT_SETUP Plots the geometry of the setup

figure;
hold on;
scatter3(scattererPoints(:, 1), scattererPoints(:, 2), scattererPoints(:, 3), 'blue', '*', 'DisplayName', 'Scatterers');
scatter3(txPoints(:, 1), txPoints(:, 2), txPoints(:, 3), 'green', 's', 'DisplayName', 'TX');
scatter3(rxPoints(:, 1), rxPoints(:, 2), rxPoints(:, 3), 'red', 'x', 'DisplayName', 'RX');
hold off;
view(45, 45);
legend;
xlabel('x (cm)');
ylabel('y (cm)');
zlabel('z (cm)');
title('Setup');
daspect([1, 1, 1]);

end

