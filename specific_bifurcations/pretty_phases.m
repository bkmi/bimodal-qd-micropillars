%% Load phases
clear;
load(strcat(data_directory(), 'specific_bifurcations/phase_collection.mat'))

%% continue
if 0
    folds2 = cell(1,numel(folds));

    for i = 1:numel(folds)
        folds2{i} = bcont(foldfuncs, folds{i}, 500, 500);
    end

    figure;
    for i = 1:numel(folds2)
        plot_branch3( folds2{i}, ...
        param, ...
        'axes_indParam', {param.feed_phase.index, param.feed_ampli.index, 'x1'}, ...
        'color', 'r', ...
        'PlotStyle', { 'LineStyle', '-', 'Marker', '.' }, ...
        'add_2_gcf', 1)
        pause(1)
    end
    
    save(strcat(datadir, 'phase_collection2'), 'folds2')
end

if exist('folds2')
    folds = folds2;
end


%% plot
figure;
for i = 1:numel(folds)
    plot_branch3( folds{i}, ...
    param, ...
    'axes_indParam', {param.feed_phase.index, param.feed_ampli.index, 'x1'}, ...
    'color', 'r', ...
    'PlotStyle', { 'LineStyle', '-', 'Marker', '.' }, ...
    'add_2_gcf', 1)
end

%% plot
figure
for i = 1:numel(phases)
    plot_branch3( phases{i}, ...
        param, ...
        'axes_indParam', {param.feed_phase.index, param.feed_ampli.index, 'x1'}, ...
        'nunst_color', phases{i}.nunst, ...
        'add_2_gcf', 1)
end