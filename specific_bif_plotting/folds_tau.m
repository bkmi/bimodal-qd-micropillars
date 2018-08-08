%% Tau seperation
%  axes = feedback amp, feedback phase, intensity
%  found in 'vfolds.m' (I think) data in 'folds_main.mat'
%  found in 'folds_low_tau.m' data in 'tau06.mat'
%  found in 'folds_low_tau.m' data in 'tau04.mat'

% Standard tau (0.8 ns)
clear;
load(strcat(specific_bif_data_dir(), ...
    'folds_main.mat'))

figure;
for i = 1:numel(folds)
    plot_branch3( folds{i}, ...
    param, ...
    'axes_indParam', {param.feed_phase.index, param.feed_ampli.index, 'x1'}, ...
    'color', 'k', ...
    'PlotStyle', { 'LineStyle', '-', 'Marker', '.' }, ...
    'add_2_gcf', 1)
end

% lower tau (0.6 ns)
clear;
load(strcat(specific_bif_data_dir(), ...
    'folds_tau06.mat'))

for i = 1:numel(folds)
    plot_branch3( folds{i}, ...
    param, ...
    'axes_indParam', {param.feed_phase.index, param.feed_ampli.index, 'x1'}, ...
    'color', 'b', ...
    'PlotStyle', { 'LineStyle', '-', 'Marker', '.' }, ...
    'add_2_gcf', 1)
end

% lowest tau (0.4 ns)
clear;
load(strcat(specific_bif_data_dir(), ...
    'folds_tau04.mat'))

for i = 1:numel(folds)
    plot_branch3( folds{i}, ...
    param, ...
    'axes_indParam', {param.feed_phase.index, param.feed_ampli.index, 'x1'}, ...
    'color', 'r', ...
    'PlotStyle', { 'LineStyle', '-', 'Marker', '.' }, ...
    'add_2_gcf', 1)
end

key = {'k', 'b', 'r'};
hold on
h = zeros(numel(key), 1);
for i = 1:numel(h)
    h(i) = plot(NaN, NaN, strcat('o', key{i}));
end
lgd = legend(h, '0.8 ns', '0.6 ns', '0.4 ns','Location','southeast');
title(lgd,'\tau_{fb}');

zlabel('Strong Field Intensity')
title('Fold Continuations')
saveas(gca, strcat(data_directory(), 'tau_fb_seperation.png'))
