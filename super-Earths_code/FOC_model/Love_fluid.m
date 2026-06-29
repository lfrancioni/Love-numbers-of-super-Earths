%-------------------------------------------------------------------------%
% --- Launching code to produce LNs and plots ----------------------------%
%-------------------------------------------------------------------------%

clear clc

%Select 'Tidal' / 'Load' / 'Shear' type of LNs
Love_type = 'Shear'; 

%Select the max harmonic degree of LNs
N_max = 25;  
n = 2:N_max;

%Select the parameter to test and the testing range 
mu_ref = 3.82e11;     % default: shear modulus of the mantle
n_steps_per_side = 3; % 3 values larger and 3 lower than the central value
variation = 0.5;      % 50% variation with respect to central value
mu_mantle_range = linspace(mu_ref * (1 - variation), ...
                           mu_ref * (1 + variation), ...
                           2 * n_steps_per_side + 1);

%Results matrixes (rows = n, columns = param. values)
h_results_inc = zeros(length(n), length(mu_mantle_range));
k_results_inc = zeros(length(n), length(mu_mantle_range));
l_results_inc = zeros(length(n), length(mu_mantle_range));

h_results_c = zeros(length(n), length(mu_mantle_range));
k_results_c = zeros(length(n), length(mu_mantle_range));
l_results_c = zeros(length(n), length(mu_mantle_range));

%LNs calculation
for j = 1:length(mu_mantle_range)
    for i = 1:length(n)
        [h, k, l] = main_incompressible(n(i), Love_type, mu_mantle_range(j));
        h_results_inc(i, j) = h;
        k_results_inc(i, j) = k;
        l_results_inc(i, j) = l;
    end
end

for j = 1:length(mu_mantle_range)
    for i = 1:length(n)
        [h, k, l] = main_compressible(n(i), Love_type, mu_mantle_range(j));
        h_results_c(i, j) = h;
        k_results_c(i, j) = k;
        l_results_c(i, j) = l;
    end
end



%PLOTTING------------------------------------------------------------------
%Colormap plots
color_style = @jet;  %Best styles: nebula - parula - jet
darken = 0.92;
cmap = darken * color_style(length(mu_mantle_range));

legend_labels = cellstr(num2str(mu_mantle_range', '\\mu = %.1e Pa'));

%Incompressible LNs
figure;
hold on;
for j = 1:length(mu_mantle_range)
    plot(n, h_results_inc(:, j), '-b','LineWidth',1,'Marker','o',...
        'MarkerSize',2, 'Color', cmap(j,:), 'MarkerFaceColor', cmap(j,:));
end
title('Love number h'''' (Incompressible)');
xlabel('Harmonic degree n');
ylabel('h_n'''' ');
grid on;
legend(legend_labels, 'Location', 'best');
ax = gca; 
ax.XScale = 'log';
xticks(2:1:10)
hold off

figure;
hold on;
for j = 1:length(mu_mantle_range)
    plot(n, k_results_inc(:, j), '-b','LineWidth',1,'Marker','o',...
        'MarkerSize',2, 'Color', cmap(j,:), 'MarkerFaceColor', cmap(j,:));
end
title('Love number k'''' (Incompressible)');
xlabel('Harmonic degree n');
ylabel('k_n'''' ');
grid on;
legend(legend_labels, 'Location', 'best');
ax = gca; 
ax.XScale = 'log';
xticks(2:1:10)
hold off

figure;
hold on;
for j = 1:length(mu_mantle_range)
    plot(n, l_results_inc(:, j), '-b','LineWidth',1,'Marker','o',...
        'MarkerSize',2, 'Color', cmap(j,:), 'MarkerFaceColor', cmap(j,:));
end
title('Love number l'''' (Incompressible)');
xlabel('Harmonic degree n');
ylabel('l_n'''' ');
grid on;
legend(legend_labels, 'Location', 'best');
ax = gca; 
ax.XScale = 'log';
xticks(2:1:10)
hold off


%Compressible LNs
figure;
hold on;
for j = 1:length(mu_mantle_range)
    plot(n, h_results_c(:, j), '-b','LineWidth',1,'Marker','o',...
        'MarkerSize',2, 'Color', cmap(j,:), 'MarkerFaceColor', cmap(j,:));
end
title('Love number h'''' (Compressible)');
xlabel('Harmonic degree n');
ylabel('h_n'''' ');
grid on;
legend(legend_labels, 'Location', 'best');
ax = gca; 
ax.XScale = 'log';
xticks(2:1:10)
hold off

figure;
hold on;
for j = 1:length(mu_mantle_range)
    plot(n, k_results_c(:, j), '-b','LineWidth',1,'Marker','o',...
        'MarkerSize',2, 'Color', cmap(j,:), 'MarkerFaceColor', cmap(j,:));
end
title('Love number k'''' (Compressible)');
xlabel('Harmonic degree n');
ylabel('k_n'''' ');
grid on;
legend(legend_labels, 'Location', 'best');
ax = gca; 
ax.XScale = 'log';
xticks(2:1:10)
hold off

figure;
hold on;
for j = 1:length(mu_mantle_range)
    plot(n, l_results_c(:, j), '-b','LineWidth',1,'Marker','o',...
        'MarkerSize',2, 'Color', cmap(j,:), 'MarkerFaceColor', cmap(j,:));
end
title('Love number l'''' (Compressible)');
xlabel('Harmonic degree n');
ylabel('l_n'''' ');
grid on;
legend(legend_labels, 'Location', 'best');
ax = gca; 
ax.XScale = 'log';
xticks(2:1:10)
hold off
