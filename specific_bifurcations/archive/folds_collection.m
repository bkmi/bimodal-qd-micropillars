%% Folds

folds = {};

% increasing intensity
folds{1} = fold_v_Weak;
folds{end+1} = highAmpFold;
folds{end+1} = foldsFromMidAmpWeakDom{3};
folds{end+1} = foldsFromMidAmpWeakDom{2};
folds{end+1} = foldsFromMidAmpWeakDom{1};
folds{end+1} = fold_v_Str;



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

%% Phases

phases = {};

% low
for i = 1:numel(phases_below_v_low)
    phases{end+1} = phases_below_v_low{i};
end

% low loops
phases{end+1} = init_phase_StrDom;
phases{end+1} = unst_phase_mid1;
phases{end+1} = init_phase_WeakDom;

% mid below double parabolas
for i = 1:numel(mid2_phases) - 1
    phases{end+1} = mid2_phases{i};
end

% mid between double parabolas
phases{end+1} = phase_high;

% high
phases{end+1} = phase_vhigh;

%% plot
figure
for i = 1:numel(phases)
    plot_branch3( phases{i}, ...
        param, ...
        'axes_indParam', {param.feed_phase.index, param.feed_ampli.index, 'x1'}, ...
        'nunst_color', phases{i}.nunst, ...
        'add_2_gcf', 1)
end
