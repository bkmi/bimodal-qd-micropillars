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
end

if exist('folds2')
    folds = folds2;
end

if 0
    save(strcat(data_directory(), 'specific_bifurcations/folds'), ...
        'foldfuncs', 'folds', 'param', 'st_bif', 'opt_inputs')
end

if 0
    % extend phase
    ok = bcont(funcs, phases{7}, 25, 0);
    ok = branch_stability(funcs, ok);
    figure
    plot_branch3( ok, ...
            param, ...
            'axes_indParam', {param.feed_phase.index, param.feed_ampli.index, 'x1'}, ...
            'nunst_color', ok.nunst, ...
            'add_2_gcf', 1)

    phases{7} = ok;

    % move 7 over by 2pi
    ph = phases{7};

    for i = 1:numel(ph.point)
        ph.point(i).parameter(param.feed_phase.index) = ph.point(i).parameter(param.feed_phase.index) - 2*pi;
    end

    figure
    for i = 1:numel(phases)
        if i == 7
            plot_branch3( ph, ...
            param, ...
            'axes_indParam', {param.feed_phase.index, param.feed_ampli.index, 'x1'}, ...
            'nunst_color', ph.nunst, ...
            'add_2_gcf', 1)
        else
            plot_branch3( phases{i}, ...
                param, ...
                'axes_indParam', {param.feed_phase.index, param.feed_ampli.index, 'x1'}, ...
                'nunst_color', phases{i}.nunst, ...
                'add_2_gcf', 1)
        end
    end

    phases{7} = ph;
end

%% extend other than 7
for i = 1:numel(phases)
    if i ~= 7
        phases{i}.parameter.min_bound = [param.feed_phase.index -18];
        phases{i}.parameter.max_bound = [param.feed_phase.index 18];
        phases{i}.parameter.max_step = [param.feed_phase.index, pi/40, ... 
                                        param.omega1.index, 0.05];
        phases{i} = bcont(funcs, phases{i}, 1000, 1000);
        phases{i} = branch_stability(funcs, phases{i});
    end
end

%% remake, other than 7
rephases = phases;

for i = 1:numel(rephases)
    if i == 6
        rephases{i}.point(end-4:end) = [];
        rephases{i}.parameter.min_bound = [param.feed_phase.index -18];
        rephases{i}.parameter.max_bound = [param.feed_phase.index 18];
        rephases{i}.parameter.max_step = [param.feed_phase.index, pi/40, ... 
                                        param.omega1.index, 0.05];
        rephases{i} = bcont(funcs, rephases{i}, 2000, 2000);
        rephases{i} = branch_stability(funcs, rephases{i});
    elseif i == 7
        rephases{i}.point(4:end) = [];
        rephases{i}.parameter.min_bound = [param.feed_phase.index -18];
        rephases{i}.parameter.max_bound = [param.feed_phase.index 18];
        rephases{i}.parameter.max_step = [param.feed_phase.index, pi/40, ... 
                                        param.omega1.index, 0.05];
        rephases{i} = bcont(funcs, rephases{i}, 200, 200);
        rephases{i} = branch_stability(funcs, rephases{i});
    elseif i == 9 || i == 10
        rephases{i}.point(4:end) = [];
        rephases{i}.parameter.min_bound = [param.feed_phase.index -18];
        rephases{i}.parameter.max_bound = [param.feed_phase.index 18];
        rephases{i}.parameter.max_step = [param.feed_phase.index, pi/40, ... 
                                        param.omega1.index, 0.05];
        rephases{i} = bcont(funcs, rephases{i}, 9000, 9000);
        rephases{i} = branch_stability(funcs, rephases{i});
    elseif i ~= 7
        rephases{i}.point(4:end) = [];
        rephases{i}.parameter.min_bound = [param.feed_phase.index -18];
        rephases{i}.parameter.max_bound = [param.feed_phase.index 18];
        rephases{i}.parameter.max_step = [param.feed_phase.index, pi/40, ... 
                                        param.omega1.index, 0.05];
        rephases{i} = bcont(funcs, rephases{i}, 5000, 5000);
        rephases{i} = branch_stability(funcs, rephases{i});
    end
end

%% save?
if 0
    save(strcat(data_directory(), 'specific_bifurcations/phases'), ...
        'funcs', 'phases', 'param', 'st_amp', 'st_phase', 'opt_inputs')
end

if 0
    save(strcat(data_directory(), 'specific_bifurcations/rephases'), ...
    'funcs', 'rephases', 'param', 'st_amp', 'st_phase', 'opt_inputs')
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
for i = 10 % numel(phases)
    plot_branch3( phases{i}, ...
        param, ...
        'axes_indParam', {param.feed_phase.index, param.feed_ampli.index, 'x1'}, ...
        'nunst_color', phases{i}.nunst, ...
        'add_2_gcf', 1)
%     pause(0.5)
end

%% plot rephaes
figure
for i = 1:numel(rephases)
    plot_branch3( rephases{i}, ...
        param, ...
        'axes_indParam', {param.feed_phase.index, param.feed_ampli.index, 'x1'}, ...
        'nunst_color', rephases{i}.nunst, ...
        'add_2_gcf', 1)
end

%% plot each
figure
for i = 1:3
    plot_branch3( rephases{i}, ...
        param, ...
        'axes_indParam', {param.feed_phase.index, param.feed_ampli.index, 'x1'}, ...
        'nunst_color', rephases{i}.nunst, ...
        'add_2_gcf', 1)
end

figure
for i = 4:6
    plot_branch3( rephases{i}, ...
        param, ...
        'axes_indParam', {param.feed_phase.index, param.feed_ampli.index, 'x1'}, ...
        'nunst_color', rephases{i}.nunst, ...
        'add_2_gcf', 1)
end

figure
for i = 7
    plot_branch3( rephases{i}, ...
        param, ...
        'axes_indParam', {param.feed_phase.index, param.feed_ampli.index, 'x1'}, ...
        'nunst_color', rephases{i}.nunst, ...
        'add_2_gcf', 1)
end

figure
for i = 8:9
    plot_branch3( rephases{i}, ...
        param, ...
        'axes_indParam', {param.feed_phase.index, param.feed_ampli.index, 'x1'}, ...
        'nunst_color', rephases{i}.nunst, ...
        'add_2_gcf', 1)
end

figure
for i = 10
    plot_branch3( rephases{i}, ...
        param, ...
        'axes_indParam', {param.feed_phase.index, param.feed_ampli.index, 'x1'}, ...
        'nunst_color', rephases{i}.nunst, ...
        'add_2_gcf', 1)
end