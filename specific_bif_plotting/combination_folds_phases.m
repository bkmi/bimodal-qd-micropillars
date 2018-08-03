%% phase continuations and folds
%  axes = feedback amp, feedback phase, intensity
%  combiniation of 'folds_main.m' and 'phase_continuations.m'

% plot folds with phases
clear;
load(strcat(specific_bif_data_dir(), ...
    'folds_main.mat'))

key = {'m', 'y', 'r', 'g', 'b', 'k'};
h = zeros(numel(folds), 1);

figure;
for i = 1:numel(folds)
    plot_branch3( folds{i}, ...
    param, ...
    'axes_indParam', {param.feed_phase.index, param.feed_ampli.index, 'x1'}, ...
    'color', key{i}, ...
    'PlotStyle', { 'LineStyle', '-', 'Marker', '.' }, ...
    'add_2_gcf', 1)
end

hold on

for i = 1:numel(h)
    h(i) = plot(NaN, NaN, strcat('o', key{i}));
end
legend(h, 'f1', 'f2', 'f3', 'f4', 'f5', 'f6');

clear;
load(strcat(specific_bif_data_dir(), ...
    'phase_continuations.mat'))

for i = 1:numel(rephases)
    plot_branch3( rephases{i}, ...
        param, ...
        'axes_indParam', {param.feed_phase.index, param.feed_ampli.index, 'x1'}, ...
        'nunst_color', rephases{i}.nunst, ...
        'add_2_gcf', 1)
end