% Demo: Reproduce towed boom profiles using elliptic integral formulas
% This script plots dimensionless profiles x/sqrt(A) vs z/sqrt(A) for
% several values of the half-included angle alpha.

clear; clc;

% Physical parameters (only the ratio T0/tau matters for shape scaling)
T0 = 1;    % tension (arbitrary units)
tau = 1;   % skin friction (same units implied)

% Half-included angles in degrees (keep below ~41 deg)
alpha_list_deg = [10, 20, 30, 40];
colors = lines(numel(alpha_list_deg));

figure('Color','w'); hold on; grid on;
for i = 1:numel(alpha_list_deg)
    alpha_deg = alpha_list_deg(i);
    prof = towed_boom_profile(alpha_deg, T0, tau, 600);

    % Plot upper branch and its mirror (lower branch)
    plot(prof.x_scaled,  prof.z_scaled,  'Color', colors(i,:), 'LineWidth', 1.8);
    plot(prof.x_scaled, -prof.z_scaled,  'Color', colors(i,:), 'LineWidth', 1.0, 'LineStyle','--');

    legends{i} = sprintf('2\\alpha = %.0f%s', 2*alpha_deg, char(176)); %#ok<SAGROW>
end
axis equal;
xlabel('x / \\surd A');
ylabel('z / \\surd A');
title('Profiles of a towed boom of logs (dimensionless scaling)');
legend(legends, 'Location','best');

% Print key scalars for the last profile as an example
prof = towed_boom_profile(alpha_list_deg(end), T0, tau, 400);
fprintf('Alpha = %.1f deg, m = %.5f\n', prof.alpha_deg, prof.m);
fprintf('A (from equilibrium) = %.6f\n', prof.A);
fprintf('A (from eq. 1.59)   = %.6f\n', prof.A_from_159);
fprintf('Total length L/sqrt(A) (eq. 1.61) = %.6f\n', prof.L_scaled);
