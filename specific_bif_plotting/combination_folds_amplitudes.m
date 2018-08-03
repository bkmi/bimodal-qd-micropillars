%% amplitude continuations and folds
%  axes = feedback amp, feedback phase, intensity
%  combiniation of 'folds_main.m' and 'amplitude_continuation_emcs.m'

%% At a fixed phase = 0
% plot folds with amplitude ecms
clear;
load(strcat(specific_bif_data_dir(), ...
    'folds_main.mat'))

wfolds = folds;
for i = 1:numel(wfolds)
    wfolds{i} = wrap_to_2pi(folds{i}, param.feed_phase.index);
end

if true
    folds = wfolds;
end

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
    'amplitude_continuation_emcs.mat'))

% plot all emcs as phase = 0
shape = size(famp_cur560_phase0);
for i = 1:shape(1)
    for j = 1:shape(2)
        if ~isempty(famp_cur560_phase0{i,j})
            plot_branch3( famp_cur560_phase0{i,j}, ...
                param, ...
                'axes_indParam', {param.feed_phase.index, param.feed_ampli.index, 'x1'}, ...
                'nunst_color', famp_cur560_phase0{i,j}.nunst, ...
                'add_2_gcf', 1)
        end
    end
end

%% At a fixed phase = pi/2
% Explain the missing folds!

% plot folds with amplitude ecms
clear;
load(strcat(specific_bif_data_dir(), ...
    'folds_main.mat'))

wfolds = folds;
for i = 1:numel(wfolds)
    wfolds{i} = wrap_to_2pi(folds{i}, param.feed_phase.index);
end

if true
    folds = wfolds;
end

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
    'amplitude_continuation_emcs.mat'))

% plot all emcs as phase = pi/2
shape = size(famp_cur560_phasePi2);
for i = 1:shape(1)
    for j = 1:shape(2)
        if ~isempty(famp_cur560_phasePi2{i,j})
            plot_branch3( famp_cur560_phasePi2{i,j}, ...
                param, ...
                'axes_indParam', {param.feed_phase.index, param.feed_ampli.index, 'x1'}, ...
                'nunst_color', famp_cur560_phasePi2{i,j}.nunst, ...
                'add_2_gcf', 1)
        end
    end
end

%% At a fixed phase = pi/4
% Explain the missing folds!

% plot folds with amplitude ecms
clear;
load(strcat(specific_bif_data_dir(), ...
    'folds_main.mat'))

wfolds = folds;
for i = 1:numel(wfolds)
    wfolds{i} = wrap_to_2pi(folds{i}, param.feed_phase.index);
end

if true
    folds = wfolds;
end

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
    'amplitude_continuation_emcs.mat'))

% plot all emcs as phase = pi/4
shape = size(famp_cur560_phasePi4);
for i = 1:shape(1)
    for j = 1:shape(2)
        if ~isempty(famp_cur560_phasePi4{i,j})
            plot_branch3( famp_cur560_phasePi4{i,j}, ...
                param, ...
                'axes_indParam', {param.feed_phase.index, param.feed_ampli.index, 'x1'}, ...
                'nunst_color', famp_cur560_phasePi4{i,j}.nunst, ...
                'add_2_gcf', 1)
        end
    end
end